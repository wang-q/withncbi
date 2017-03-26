#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Template;
use Path::Tiny;
use List::MoreUtils qw(uniq zip);
use Text::CSV_XS;

use Bio::DB::Taxonomy;

use AlignDB::Stopwatch;

use lib "$FindBin::Bin/../lib";
use MyUtil qw(wgs_worker);

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

=head1 NAME

wgs_prep.pl - prepare WGS materials

=head1 SYNOPSIS

    perl wgs_prep.pl [options]
      Options:
        --help      -?          brief help message

        -i, -f, --file      tab seperated file containing wgs prefix and name
        -o, -d, --dir       output dir
        --aria2     -a          url file is for aria2
        --fix                   sometimes WGS records miss assigning strain id
        --nofix         @STR    Skip some strains
        --csvonly               only generate the csv file

    perl wgs_prep.pl -a -f trichoderma.tsv -o WGS
    Three files will be generated.
    trichoderma.csv
    trichoderma.url.txt
    trichoderma.data.txt

=cut

my $td_dir = path( $Config->{path}{td} )->stringify;    # taxdmp

# for unrecorded strains, give them arbitrary ids
my $arbitrary = 100_000_000;

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'i|f|file=s' => \my $file_input,
    'o|d|dir=s'  => \my $dir_output,
    'a|aria2'    => \my $aria2,
    'fix'        => \my $fix_strain,
    'nofix=s'    => \my @nofix,
    'csvonly'    => \my $csvonly,
) or HelpMessage(1);

$dir_output = "." unless $dir_output;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Prepare NCBI WGS");

my $taxon_db = Bio::DB::Taxonomy->new( -source => 'entrez', );

#----------------------------#
# Read
#----------------------------#
$stopwatch->block_message("Load $file_input.");
my $basename = path($file_input)->basename( ".txt", ".tab", ".tsv" );

my $wgsid_of = {};
my @orig_orders;
{
    my @lines = path($file_input)->lines;
    for my $line (@lines) {
        chomp $line;
        $line =~ /^#/ and next;
        my ( $name, $prefix ) = split /\t/, $line;
        $prefix or next;
        $wgsid_of->{$name} = $prefix;
        push @orig_orders, $name;
    }
}

#----------------------------#
# scraper
#----------------------------#
$stopwatch->block_message("Scrapping NCBI WGS...");
my $master = {};
{
    for my $key (@orig_orders) {
        print "$key\n";
        my $prefix = $wgsid_of->{$key};
        my $info   = wgs_worker($prefix);
        $info->{name} = $key;
        $master->{$key} = $info;
    }

    print "\n", "=" x 20, "\n";
    print "Finish scrapping\n";
}

#----------------------------#
# csv and url
#----------------------------#
$stopwatch->block_message("Generate .csv for info and .url.txt for downloading ");
{
    mkdir $dir_output unless -d $dir_output;

    my $csv = Text::CSV_XS->new( { binary => 1 } )
        or die "Cannot use CSV: " . Text::CSV_XS->error_diag;
    $csv->eol("\n");

    my $file_csv = path( $dir_output, "$basename.csv" )->stringify;

    open my $csv_fh, ">", $file_csv;

    my @columns = (
        'prefix',                'taxon_id',     'name',         'Organism',
        'Biosource',             'BioProject',   'Keywords',     'Genome_Coverage',
        'Sequencing_Technology', '#_of_Contigs', 'Total_length', 'Assembly_Method',
        'Assembly_Name',         'Update_date',  'pubmed',
    );

    $csv->print( $csv_fh, \@columns );

    for my $key (@orig_orders) {

        # Don't use hashref here, because I want use hash slices.
        my %info = %{ $master->{$key} };

        if ( !$csvonly and $fix_strain ) {
            if ( grep { $_ eq $key } @nofix ) {
                print "Skip $info{name} as you don't want fix it\n";
            }
            elsif ( $info{Organism} =~ /$info{Biosource}/ ) {
                print "Don't need fixing for $info{name}\n";
            }
            else {
                # Sometimes the uploader didn't create a new strain, assign its own id and marked
                # species name as strain name.
                # So try looking up this strain in taxonomy dumps
                print "Fix strain taxon info for $info{name}\n";
                $info{Organism} = $info{Organism} . " " . $info{Biosource};
                my $node = $taxon_db->get_taxon( -name => $info{Organism} );
                if ( !$node ) {
                    print " " x 4, "Can't find taxon for $info{name}\n";
                    $arbitrary++;
                    print " " x 4, "Give it arbitrary id as $arbitrary\n";
                    $info{original_id} = $info{taxon_id};
                    $info{taxon_id}    = $arbitrary;
                }
                else {
                    $info{taxon_id} = $node->id;
                }
            }
        }
        $csv->print( $csv_fh, [ @info{@columns} ] );
        $master->{$key} = \%info;

    }
    close $csv_fh;

    print "\n", "=" x 20, "\n";
    print ".csv generated.\n";

    if ( !$csvonly ) {
        my $file_url = path( $dir_output, "$basename.url.txt" )->stringify;
        open my $url_fh, ">", $file_url;
        for my $key (@orig_orders) {

            my %info = %{ $master->{$key} };

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
        close $url_fh;
        if ($aria2) {
            print "# Run the follow cmd to download with aria2c\n";
            print "aria2c -UWget -x 6 -s 3 -c -i $file_url\n";
            print "\n";
            print "# Use the following cmd to check .gz files\n";
            print "find $dir_output -name \"*.gz\" | xargs gzip -t \n";
        }
        else {
            print "Download urls in file $file_url\n";
        }
    }
}

#----------------------------#
# @data yaml
#----------------------------#
if ( !$csvonly ) {
    $stopwatch->block_message("Generate .data.yml");

    my $file_data = path( $dir_output, "$basename.data.yml" )->stringify;

    my $text = <<'EOF';
---
data:
[% FOREACH name IN names -%]
  - taxon: [% master.$name.taxon_id %]
    name: [% master.$name.name %]
    sciname: [% master.$name.Organism %]
    prefix: [% master.$name.prefix %]
    coverage: [% master.$name.Genome_Coverage %] [% master.$name.Sequencing_Technology %]
[% IF master.$name.original_id -%]
    original_id: [% master.$name.original_id %]
[% END -%]
[% END -%]

EOF
    my $tt = Template->new;
    $tt->process( \$text, { names => [@orig_orders], master => $master, }, $file_data )
        or die Template->error;

    print ".data.txt generated.\n";
}

$stopwatch->end_message;

exit;

__END__
