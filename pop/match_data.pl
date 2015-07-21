#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use File::Find::Rule;
use File::Spec;
use String::Compare;

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home wgs_worker);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

my $file_input;
my $file_output;
my $dir_scan;

my $match_field = 'name';

my @name_rules = ( "*.fa", "*.fa.gz" );
my $pattern;

my %skip;
my @per_seq;

my %other_opts;

my $yes;

my $man  = 0;
my $help = 0;

GetOptions(
    'help'        => \$help,
    'man'         => \$man,
    'i|input=s'   => \$file_input,
    'o|output=s'  => \$file_output,
    'd|dir=s'     => \$dir_scan,
    'm|match=s'   => \$match_field,
    'r|rule=s'    => \@name_rules,
    'p|pattern=s' => \$pattern,
    'opt=s'       => \%other_opts,
    'skip=s'      => \%skip,
    'per_seq=s'   => \@per_seq,
    'y|yes'       => \$yes,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ($file_output) {
    $file_output = $file_input;
}

my %per_seq = map { $_ => 1 } @per_seq;

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
$stopwatch->start_message;

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
    for my $item (@data) {
        printf "Matching [%s]\n", $item->{name};

        # match the most similar name
        my ($fasta) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( basename($_), $item->{$match_field} ) ] }
            @files;
        $item->{fasta} = $fasta;
        printf " " x 4 . "%s => %s\n", $item->{$match_field}, $item->{fasta};
    }
}

$stopwatch->block_message("Mark flags");
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

# Move skipped entries to tail
my @data_sort = grep { !exists $_->{skip} } @data;
push @data_sort, grep { exists $_->{skip} } @data;

# Move per_seq entries to head
my @data_sort2 = grep { exists $_->{per_seq} } @data_sort;
push @data_sort2, grep { !exists $_->{per_seq} } @data_sort;
$yml->{data} = \@data_sort2;

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
            print
                "[$file_output] exists and may contain manual added infomaton. Don't overwrite it.\n";
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

=head1 NAME

match_data.pl - find matched files for each @data entry in YAML and store extra options.

=head1 SYNOPSIS

    perl match_data.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -i, --input         input yaml
        -o, --output        output yaml
        -d, --dir           where sequence files live
        -m, --match         key name of each @data entry. or prefix, sciname...
        -r, --rule          @, File::Find::Rule, '*.fsa_nt.gz' for NCBI WGS
        -p, --pattern       For ensembl, 'dna.toplevel'
        --opt               %, Other options for running pop
        --skip              %, skip this strain
        --per_seq           @, split fasta by names, target or good assembles
        -y, --yes           Overwrite existing YAML file

=cut
