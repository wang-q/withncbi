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

my %other_opts;

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
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ($file_output) {
    $file_output = $file_input;
}

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
$stopwatch->start_message("Load YAML");

my $yml  = LoadFile($file_input);
my @data = @{ $yml->{data} };

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
        map { [ $_, compare( basename($_), $item->{$match_field} ) ] } @files;
    $item->{fasta} = $fasta;
    printf " " x 4 . "%s => %s\n", $item->{$match_field}, $item->{fasta};
}

# Move skipped entries to tail
my @data_sort = grep { !exists $_->{skip} } @data;
push @data_sort, grep { exists $_->{skip} } @data;
$yml->{data} = \@data_sort;

$stopwatch->block_message("Other options");
for my $key ( sort keys %other_opts ) {
    $yml->{$key} = $other_opts{$key};
}

$stopwatch->block_message("Write YAML [$file_output]");
DumpFile( $file_output, $yml );

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

    perl match_data.pl \
        -i ~/data/alignment/trichoderma/WGS/trichoderma.data.yml \
        -o trichoderma_test.yml \
        -d ~/data/alignment/trichoderma/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=trichoderma \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/trichoderma' \
        --opt pl_dir='~/Scripts' \
        --opt parallel=4

=cut
