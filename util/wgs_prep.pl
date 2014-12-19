#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use Template;
use List::MoreUtils qw(uniq zip);
use Text::CSV_XS;
use File::Basename;
use File::Slurp;
use File::Spec;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(wgs_worker);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $file_input;
my $dir_output;
my $aria2;    # generate a aria2 input file

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'i|f|file=s' => \$file_input,
    'o|d|dir=s'  => \$dir_output,
    'a|aria2'    => \$aria2,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$dir_output = "." unless $dir_output;

#----------------------------------------------------------#
# Read
#----------------------------------------------------------#
my $basename = basename( $file_input, ".txt", ".tab", ".tsv" );

my $id_of = {};
{
    my @lines = read_file($file_input);
    for my $line (@lines) {
        $line =~ /^#/ and next;
        my ( $name, $prefix ) = split /\t/, $line;
        $prefix or next;
        $id_of->{$name} = $prefix;
    }
}

#----------------------------------------------------------#
# scraper
#----------------------------------------------------------#
my $master = {};
{
    for my $key ( sort keys %{$id_of} ) {
        print "$key\n";
        my $prefix = $id_of->{$key};
        my $info   = wgs_worker($prefix);
        $info->{name} = $key;
        $master->{$key} = $info;
    }

    print "\n", "=" x 20, "\n";
    print "Finish scrapping\n";
}

#----------------------------------------------------------#
# csv and url
#----------------------------------------------------------#
{
    mkdir $dir_output unless -d $dir_output;

    my $csv = Text::CSV_XS->new( { binary => 1 } )
        or die "Cannot use CSV: " . Text::CSV_XS->error_diag;
    $csv->eol("\n");

    my $file_csv = File::Spec->catfile( $dir_output, "$basename.csv" );
    my $file_url = File::Spec->catfile( $dir_output, "$basename.url.txt" );

    open my $csv_fh, ">", $file_csv;
    open my $url_fh, ">", $file_url;

    my @columns = (
        'prefix',          'taxon_id',
        'name',            'Organism',
        'BioProject',      'Keywords',
        'Genome_Coverage', 'Sequencing_Technology',
        '#_of_Contigs',    'Total_length',
        'Assembly_Method', 'Assembly_Name',
        'Update_date',     'pubmed',
    );

    $csv->print( $csv_fh, \@columns );

    for my $key ( sort keys %{$id_of} ) {
        my %info = %{ $master->{$key} };
        my $row  = [ @info{@columns} ];
        $csv->print( $csv_fh, $row );

        my @downloads = @{ $info{download} };
        if ($aria2) {
            for (@downloads) {
                print {$url_fh} $_, "\n";
                print {$url_fh} "  dir=", "$dir_output/$key", "\n";
            }
        }
        else {
            print {$url_fh} $_, "\n" for @downloads;
        }
    }

    close $csv_fh;
    close $url_fh;

    print "\n", "=" x 20, "\n";
    print ".csv generated.\n";

    print "\n", "=" x 20, "\n";
    if ($aria2) {
        print "# Run the follow cmd to download with aria2c\n";
        print "aria2c -x 6 -s 3 -c -i $file_url\n";
    }
    else {
        print "Download files in $file_url\n";
    }
    print "# Use the following cmd to check .gz files\n";
    print "find $dir_output -name \"*.gz\" | xargs gzip -t \n";
}

#----------------------------------------------------------#
# @data template
#----------------------------------------------------------#
{
    my $file_data = File::Spec->catfile( $dir_output, "$basename.data.txt" );

    my $text = <<'EOF';
my @data = (
[% FOREACH name IN names -%]
    {   taxon    => [% master.$name.taxon_id %],
        name     => "[% master.$name.name %]",
        sciname  => "[% master.$name.Organism %]",
        prefix   => "[% master.$name.prefix %]",
        coverage => "[% master.$name.Genome_Coverage %] [% master.$name.Sequencing_Technology %]",
    },
[% END -%]
);
EOF
    my $tt = Template->new;
    $tt->process( \$text,
        { names => [ sort keys %{$id_of} ], master => $master, }, $file_data )
        or die Template->error;

    print "\n", "=" x 20, "\n";
    print ".data.txt generated.\n";
}

exit;

=head1 NAME

    wgs_prep.pl - prepare for wgs

=head1 SYNOPSIS

    perl wgs_prep.pl -a -f trichoderma.tsv -o WGS

    wgs_prep.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -i, -f, --file      tab seperated file containing wgs prefix and name
        -o, -d, --dir       output dir

    Three files will be generated.
    trichoderma.csv
    trichoderma.url.txt
    trichoderma.data.txt

=cut
