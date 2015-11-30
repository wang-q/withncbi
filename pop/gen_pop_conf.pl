#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use File::Find::Rule;
use File::Spec;
use String::Compare;
use List::MoreUtils qw(zip);
use Set::Scalar;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

gen_pop_conf.pl - for each @data entries in YAML, find matched files, check parameters and store other options.

=head1 SYNOPSIS

    perl gen_pop_conf.pl <-i data.yml> [options]
      Options:
        --help      -?          brief help message
        --input     -i  STR     input yaml
        --output    -o  STR     output yaml
        --dir       -d  STR     Where sequence files live
        --match     -m  STR     Key of each @data entry. name, prefix or sciname...
        --rule      -r  @STR    File::Find::Rule, '*.fsa_nt.gz' for NCBI WGS
        --pattern   -p  STR     For ensembl, 'dna.toplevel'
        --opt       STR=STR     Other options for running pop
        --skip      STR=STR     Skip this strain
        --per_seq       @STR    Split fasta by names, target or good assembles
        --arbitrary     INT     For unrecorded strains, give them arbitrary ids, default is [100_000_000]
        --dd            STR     Where downloaded files live
                                Default is ".".
        --download      @STR    Add entries to @data which were download previously
                                'name=Scer_S288c;taxon=559292;sciname=Saccharomyces cerevisiae S288c'
        --plan          @STR    Add alignment plans
                                'name=four_way;t=Scer_S288c;qs=Sbou_ATCC_MYA_796,Spar_NRRL_Y_17217,Spas_CBS_1483'
        --yes       -y          Overwrite existing YAML file

=cut

#

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'input|i=s'  => \my $file_input,
    'output|o=s' => \my $file_output,
    'dir|d=s'    => \my $dir_scan,
    'match|m=s' => \( my $match_field = 'name' ),
    'rule|r=s' => \my @name_rules,
    'pattern|p=s' => \my $pattern,
    'opt=s'       => \my %other_opts,
    'skip=s'      => \my %skip,
    'per_seq=s'   => \my @per_seq,
    'arbitrary=i' => \( my $arbitrary = 100_000_000 ),
    'dd=s'        => \( my $dir_download = "." ),
    'download=s'  => \my @downloaded,
    'plan=s'      => \my @plan,
    'yes|y'       => \my $yes,
) or HelpMessage(1);

die "Need a YAML file" unless $file_input;

unless ($file_output) {
    $file_output = $file_input;
}

my %per_seq = map { $_ => 1 } @per_seq;

unless ( scalar @name_rules ) {
    @name_rules = ( "*.fa", "*.fa.gz" );
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
$stopwatch->start_message;

$stopwatch->block_message("Load YAML");
my $yml  = LoadFile($file_input);
my @data = @{ $yml->{data} };

if ( defined $dir_scan and -d $dir_scan ) {
    $stopwatch->block_message("Scan $dir_scan");
    my @files = File::Find::Rule->file->name(@name_rules)->in($dir_scan);
    if ($pattern) {
        my $rx = quotemeta $pattern;
        @files = grep {m{$rx}} @files;
    }

    $stopwatch->block_message("Match field");

    my @name_rules_copy = map { s/\*//; $_ } @name_rules;
    my @strips = map { basename( $_, @name_rules_copy ) } @files;
    my $file_of = { zip( @strips, @files ) };

    for my $item (@data) {
        printf "Matching [%s]\n", $item->{name};

        # match the most similar name
        my ($fasta) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( $_, $item->{$match_field} ) ] }
            keys %{$file_of};
        $item->{fasta} = $file_of->{$fasta};
        printf " " x 4 . "%s => %s => %s\n", $item->{$match_field}, $fasta, $item->{fasta};

        if ( index( lc( $item->{fasta} ), lc( $item->{name} ) ) == -1 ) {
            printf " " x 4 . "[%s] with [%s] matches to [%s]\n", $item->{name},
                $item->{$match_field}, $item->{fasta};
            print " " x 4 . "Filtering with prefix and try again.\n";
            die "Match errors. Please check.\n";
        }
    }
}

$stopwatch->block_message("Mark flags");
my $name_set = Set::Scalar->new;
$name_set->insert( $_->{name} ) for @data;
for my $name ( keys %per_seq, keys %skip ) {
    if ( !$name_set->has($name) ) {
        die "Check you --skip or --per_seq for [$name], which isn't present in YAML-data-names.\n";
    }
}

for my $item (@data) {
    if ( $skip{ $item->{name} } ) {
        printf "[%s] Mark flag 'skip'\n", $item->{name};
        $item->{skip} = $skip{ $item->{name} };
    }
    if ( $per_seq{ $item->{name} } ) {
        printf "[%s] Mark flag 'per_seq'\n", $item->{name};
        $item->{per_seq} = 1;
    }
}

# arbitrary ids
for my $item (@data) {
    if ( !exists $item->{taxon} ) {
        $arbitrary++;
        $item->{taxon} = $arbitrary;
    }
}

$stopwatch->block_message('Sort @data');

# Move skipped entries to tail
my @data_sort = grep { !exists $_->{skip} } @data;
push @data_sort, grep { exists $_->{skip} } @data;

# Move per_seq entries to head
my @data_sort2 = grep { exists $_->{per_seq} } @data_sort;
push @data_sort2, grep { !exists $_->{per_seq} } @data_sort;

# 'name=Scer_S288c,taxon=559292,sciname=Saccharomyces cerevisiae S288c'
if ( scalar @downloaded ) {
    $stopwatch->block_message('Pre downloaded');

    for my $entry ( reverse @downloaded ) {
        my %hash = map { split /=/ } ( split /;/, $entry );
        printf "Inject downloaded %s\n", $hash{name};

        $hash{downloaded} = 1;
        my $dir = File::Spec->catdir( $dir_download, $hash{name} );
        if ( -d $dir ) {
            $hash{pre_dir} = $dir;
            printf "%s => %s\n", $hash{name}, $hash{pre_dir};
        }
        else {
            printf
                "%s doesn't exist, make sure %s/ exists in working directory.\n",
                $dir, $hash{name};
        }

        if ( $name_set->has( $hash{name} ) ) {

            # replace existing one
            @data_sort2 = grep { $_->{name} ne $hash{name} } @data_sort2;
        }
        else {
            # add to $name_set
            $name_set->insert( $hash{name} );
        }

        # to the head
        unshift @data_sort2, \%hash;
    }
}
$yml->{data} = \@data_sort2;

if ( scalar @plan ) {
    $stopwatch->block_message("Alignment plans");
    my @ary;
    for my $entry (@plan) {
        my %hash = map { split /=/ } ( split /;/, $entry );
        $hash{qs} = [ split /,/, $hash{qs} ];

        for my $name ( $hash{t}, $hash{o}, @{ $hash{qs} } ) {
            next unless defined $name;
            if ( !$name_set->has($name) ) {
                printf "In plan [%s]\n", $hash{name};
                die "Please check for [$name], which isn't present in YAML-data-names.\n";
            }
        }

        printf "Inject plan %s\n", $hash{name};
        push @ary, \%hash;
    }
    $yml->{plans} = \@ary;
}

$stopwatch->block_message("Other options");
for my $key ( sort keys %other_opts ) {
    $yml->{$key} = $other_opts{$key};
}

if ( -e $file_output ) {
    if ( $file_input ne $file_output ) {
        if ($yes) {
            $stopwatch->block_message("Write YAML [$file_output]");
            DumpFile( $file_output, $yml );
        }
        else {
            $stopwatch->block_message("NOTICE");
            print
                "[$file_output] exists and may contain manually added infomaton. Don't overwrite it.\n";
        }
    }
    else {
        $stopwatch->block_message("Update YAML [$file_output]");
        DumpFile( $file_output, $yml );
    }
}
else {
    $stopwatch->block_message("Write YAML [$file_output]");
    DumpFile( $file_output, $yml );
}

$stopwatch->end_message;

exit;

__END__
