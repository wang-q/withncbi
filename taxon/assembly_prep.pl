#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw();
use YAML::Syck qw();
use Path::Tiny qw();

use Template;
use Text::CSV_XS;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();

=head1 NAME

assembly_prep.pl - prepare ASSEMBLY materials

=head1 SYNOPSIS

    perl assembly_prep.pl [options]
      Options:
        --help      -?          brief help message

        --file, -f      STR     tab seperated file containing wgs prefix and name
        --outdir, -o    STR     output dir
        --csvonly               only generate the csv file

    perl assembly_prep.pl -f trichoderma.assembly.tsv -o ASSEMBLY

    #name   ftp_path    organism    assembly_level

    Three files will be generated.
    trichoderma.csv
    trichoderma.url.txt
    trichoderma.data.txt

=cut

Getopt::Long::GetOptions(
    'help|?'     => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'   => \my $infile,
    'outdir|o=s' => \( my $outdir = "." ),
    'csvonly'    => \my $csvonly,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Prepare NCBI ASSEMBLY");

#----------------------------#
# Read
#----------------------------#
$stopwatch->block_message("Load $infile.");
my $basename = Path::Tiny::path($infile)->basename( ".txt", ".tab", ".tsv" );

my $ftp_of = {};
my @orig_orders;
{
    my @lines = Path::Tiny::path($infile)->lines;
    for my $line (@lines) {
        chomp $line;
        $line =~ /^#/ and next;
        my ( $name, $ftp ) = split /\t/, $line;
        $ftp or next;
        $ftp_of->{$name} = $ftp;
        push @orig_orders, $name;
    }
}

#----------------------------#
# rsync .sh
#----------------------------#
$stopwatch->block_message("Generate rsync .sh for downloading ");

Path::Tiny::path($outdir)->mkpath();
my $file_sh = Path::Tiny::path( $outdir, "$basename.sh" );
$file_sh->remove if $file_sh->is_file;

$file_sh->append(
    <<'EOF'
#!/bin/bash

BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${BASE_DIR}

EOF
);

for my $key (@orig_orders) {
    my $ftp = $ftp_of->{$key};

    # ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/167/675/GCA_000167675.2_v2.0
    # ftp.ncbi.nlm.nih.gov::genomes/all/GCA/000/167/675/GCA_000167675.2_v2.0

    my $rsync = $ftp;
    $rsync =~ s/ftp:\/\/ftp.ncbi.nlm.nih.gov\//ftp.ncbi.nlm.nih.gov::/;
    if ( $rsync eq $ftp ) {
        die "Check the ftp url: $ftp\n";
    }

    $file_sh->append(
        <<"EOF"
echo >&2
echo >&2 "==> $key"
mkdir -p $key
rsync -avP $rsync/ $key/

EOF
    );
}

#{
#    mkdir $outdir unless -d $outdir;
#
#    my $csv = Text::CSV_XS->new( { binary => 1 } )
#        or die "Cannot use CSV: " . Text::CSV_XS->error_diag;
#    $csv->eol("\n");
#
#    my $file_csv = Path::Tiny::path( $outdir, "$basename.csv" )->stringify;
#
#    open my $csv_fh, ">", $file_csv;
#
#    my @columns = (
#        'prefix',                'taxon_id',     'name',         'Organism',
#        'Biosource',             'BioProject',   'Keywords',     'Genome_Coverage',
#        'Sequencing_Technology', '#_of_Contigs', 'Total_length', 'Assembly_Method',
#        'Assembly_Name',         'Update_date',  'pubmed',
#    );
#
#    $csv->print( $csv_fh, \@columns );
#
#    for my $key (@orig_orders) {
#
#        # Don't use hashref here, because I want use hash slices.
#        my %info = %{ $master->{$key} };
#
#        if ( !$csvonly and $fix_strain ) {
#            if ( grep { $_ eq $key } @nofix ) {
#                print "Skip $info{name} as you don't want fix it\n";
#            }
#            elsif ( $info{Organism} =~ /$info{Biosource}/ ) {
#                print "Don't need fixing for $info{name}\n";
#            }
#            else {
#                # Sometimes the uploader didn't create a new strain, assign its own id and marked
#                # species name as strain name.
#                # So try looking up this strain in taxonomy dumps
#                print "Fix strain taxon info for $info{name}\n";
#                $info{Organism} = $info{Organism} . " " . $info{Biosource};
#
#                my $node;
#                eval { $node = $taxon_db->get_taxon( -name => $info{Organism} ); };
#                if ( $@ or !$node ) {
#                    print " " x 4, "Can't find taxon for $info{name}\n";
#                    $arbitrary++;
#                    print " " x 4, "Give it arbitrary id as $arbitrary\n";
#                    $info{original_id} = $info{taxon_id};
#                    $info{taxon_id}    = $arbitrary;
#                }
#                else {
#                    $info{taxon_id} = $node->id;
#                }
#            }
#        }
#        $csv->print( $csv_fh, [ @info{@columns} ] );
#        $master->{$key} = \%info;
#
#    }
#    close $csv_fh;
#
#    print "\n", "=" x 20, "\n";
#    print ".csv generated.\n";
#
#    if ( !$csvonly ) {
#        my $file_url = Path::Tiny::path( $outdir, "$basename.url.txt" )->stringify;
#        open my $url_fh, ">", $file_url;
#        for my $key (@orig_orders) {
#
#            my %info = %{ $master->{$key} };
#
#            my @downloads = @{ $info{download} };
#            if ($aria2) {
#                for (@downloads) {
#                    print {$url_fh} $_, "\n";
#                    print {$url_fh} "  dir=", "$outdir/$key", "\n";
#                }
#            }
#            else {
#                print {$url_fh} $_, "\n" for @downloads;
#            }
#
#        }
#        close $url_fh;
#        if ($aria2) {
#            print "# Run the follow cmd to download with aria2c\n";
#            print "aria2c -UWget -x 6 -s 3 -c -i $file_url\n";
#            print "\n";
#            print "# Use the following cmd to check .gz files\n";
#            print "find $outdir -name \"*.gz\" | xargs gzip -t \n";
#        }
#        else {
#            print "Downloading urls in file $file_url\n";
#        }
#    }
#}

#----------------------------#
# @data yaml
#----------------------------#
#if ( !$csvonly ) {
#    $stopwatch->block_message("Generate .data.yml");
#
#    my $file_data = Path::Tiny::path( $outdir, "$basename.data.yml" )->stringify;
#
#    my $text = <<'EOF';
#---
#data:
#[% FOREACH name IN names -%]
#  - taxon: [% master.$name.taxon_id %]
#    name: [% master.$name.name %]
#    sciname: [% master.$name.Organism %]
#    prefix: [% master.$name.prefix %]
#    coverage: [% master.$name.Genome_Coverage %] [% master.$name.Sequencing_Technology %]
#[% IF master.$name.original_id -%]
#    original_id: [% master.$name.original_id %]
#[% END -%]
#[% END -%]
#
#EOF
#    my $tt = Template->new;
#    $tt->process( \$text, { names => [@orig_orders], master => $master, }, $file_data )
#        or die Template->error;
#
#    print ".data.yml generated.\n";
#}

$stopwatch->end_message;

exit;

__END__
