#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use Path::Tiny;
use Text::CSV_XS;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

batch_get_seq.pl - retrieve all sequences listed in a file

=head1 SYNOPSIS

    perl batch_get_seq.pl <-f seq.csv> [options]
      Options:
        --help      -?          brief help message
        --file      -f  STR     input csv file
        --local     -l  STR     find files from this local directory first
        --pure      -p          remove fasta sequences from generated gff files

=head1 CSV file format

Example:

    strain_name,accession,strain_taxon_id,seq_name # This line is needed
    S288c,NC_001133,559292,I
    S288c,NC_001134,559292,II
    S288c,NC_001135,559292,III

Columns:

    strain_name     =>  files will be stored in a directory named after this
    accession       =>  NCBI sequence accession
    strain_taxon_id =>  optional, not needed
    seq_name        =>  optional, with -r, rename filename and sequence header

=head1 EXAMPLE

    cd ~/data/alignment/yeast_genome
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl -f yeast_name_seq.csv 2>&1 |
        tee yeast_name_seq.log

=cut

GetOptions(
    'help|?'    => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s'  => \my $in_file,
    'local|l=s' => \my $local,
) or Getopt::Long::HelpMessage(1);

die "Provide a .csv file\n" unless defined $in_file and -e $in_file;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my ( @fna_files, @gff_files );
if ($local) {
    $local = path($local)->absolute->stringify;
    print "Reading file list from [$local]\n";
    @fna_files = File::Find::Rule->file->name('*.fna')->in($local);
    @gff_files = File::Find::Rule->file->name('*.gff')->in($local);
    printf "Get [%d] fna files\n", scalar @fna_files;
}

my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
open my $csv_fh, "<", $in_file or die "$in_file: $!";
$csv->getline($csv_fh);    # bypass title line

#----------------------------------------------------------#
# Download
#----------------------------------------------------------#
while ( my $row = $csv->getline($csv_fh) ) {
    my $id  = $row->[0];
    my $acc = $row->[1];

    # replace non-alphanumeric chars
    $id =~ s/[\W]+/_/g;

    print "==> id: [$id]\tseq: [$acc]\n";
    if ( -e "$id/$acc.gff" and -e "$id/$acc.fasta" ) {
        print "Sequence [$id/$acc] exists, next\n\n";
        next;
    }

    if ( !$local ) {
        print "Fetch sequences from NCBI\n";
        system "perl $FindBin::RealBin/get_seq.pl $acc $id";
    }
    else {
        print "Try finding sequences from local disk\n";
        my $id_dir = path($id);
        if ( !-e $id_dir ) {
            $id_dir->mkpath;
        }

        my ($fna_file) = grep {/$acc/} @fna_files;
        if ( !defined $fna_file ) {
            print "Fetch sequences from NCBI\n";
            system "perl $FindBin::RealBin/get_seq.pl $acc $id";
        }
        else {
            path($fna_file)->copy("$id/$acc.fasta");

            my ($gff_file) = grep {/$acc/} @gff_files;
            path($gff_file)->copy("$id/$acc.gff");
        }
    }

    # replace contents from .gb and .fasta
    if ( defined $row->[3] ) {
        my $seq_name = $row->[3];
        if ( -e "$id/$acc.gb" ) {
            system "perl -i -nlp -e '/^(?:LOCUS)/ and s/$acc/$seq_name/' $id/$acc.gb";
        }
        if ( -e "$id/$acc.fasta" ) {
            system "perl -i -nlp -e '/^(?:>)/ and s/.+/>$seq_name/' $id/$acc.fasta";
        }
    }

    if ( -e "$id/$acc.gb" ) {
        system "perl $FindBin::Bin/bp_genbank2gff3.pl $id/$acc.gb";
        Path::Tiny::path("$acc.gb.gff")->move("$id/$acc.gff");
        system "perl -i -nlp -e '/^\#\#FASTA/ and last' $id/$acc.gff";
    }

    if ( defined $row->[3] ) {
        my $seq_name = $row->[3];
        print " " x 4 . "Rename .fasta and .gff\n";
        path("$id/$acc.fasta")->move("$id/$seq_name.fasta");
        path("$id/$acc.gff")->move("$id/$seq_name.gff");
    }

    print "\n\n";
}

close $csv_fh;

__END__
