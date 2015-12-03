#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Text::CSV_XS;
use File::Copy;

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
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl -p -f yeast_name_seq.csv 2>&1 | tee yeast_name_seq.log

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \my $in_file,
    'pure|p'   => \my $pure_gff,
) or HelpMessage(1);

die "Provide a .csv file\n" unless defined $in_file and -e $in_file;

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#
my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
open my $csv_fh, "<", $in_file or die "$in_file: $!";
$csv->getline($csv_fh);    # bypass title line

#----------------------------------------------------------#
# Download
#----------------------------------------------------------#
while ( my $row = $csv->getline($csv_fh) ) {
    my $id  = $row->[0];
    my $seq = $row->[1];

    # replace non-alphanumeric chars
    $id =~ s/[\W]+/_/g;

    print "id: [$id]\tseq: [$seq]\n";
    if ( -e "$id/$seq.gb" ) {
        print "Sequence [$id/$seq.gb] exists, next\n";
        next;
    }

    # download
    system "perl $FindBin::Bin/get_seq.pl $seq $id";

    # replace contents from .gb and .fasta
    if ( defined $row->[3] ) {
        my $seq_name = $row->[3];
        print " " x 4 . "Replace locus $seq to $seq_name\n";

        system "perl -i -nlp -e '/^(?:LOCUS)/ and s/$seq/$seq_name/' $id/$seq.gb";
        system "perl -i -nlp -e '/^(?:>)/ and s/.+/>$seq_name/' $id/$seq.fasta";
    }
    else {
        print " " x 4 . "Keep original names.\n";
    }

    system "perl $FindBin::Bin/bp_genbank2gff3.pl $id/$seq.gb";

    move( "$seq.gb.gff", "$id/$seq.gff" );
    if ($pure_gff) {
        system "perl -i -nlp -e '/^\#\#FASTA/ and last' $id/$seq.gff";
    }

    if ( defined $row->[3] ) {
        my $seq_name = $row->[3];
        print " " x 4 . "Rename .fasta and .gff\n";
        move( "$id/$seq.fasta", "$id/$seq_name.fasta" );
        move( "$id/$seq.gff",   "$id/$seq_name.gff" );
    }

    print "\n\n";
}

close $csv_fh;

__END__
