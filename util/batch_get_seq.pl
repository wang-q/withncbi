#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;

use Text::CSV_XS;
use File::Copy;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_file;

my $rename;    # rename source to seq name
my $pure_gff; # remove fasta sequences from generated gff files

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'   => \$help,
    'man'      => \$man,
    'f|file=s' => \$in_file,
    'r|rename' => \$rename,
    'p|pure' => \$pure_gff,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

die "Provide a csv file\n" unless defined $in_file and -e $in_file;

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

    print "id: [$id]\tseq: [$seq]\n";
    if ( -e "$id/$seq.gb" ) {
        print "Sequence [$id/$seq.gb] exists, next\n";
        next;
    }

    # download
    system "perl $FindBin::Bin/../util/get_seq.pl $seq $id";

    # replace contents from .gb and .fasta
    if ( $rename and defined $row->[3] ) {
        my $seq_name = $row->[3];
        print " " x 4 . "Replace locus $seq to $seq_name\n";

        system
            "perl -i -nlp -e '/^(?:LOCUS)/ and s/$seq/$seq_name/' $id/$seq.gb";
        system "perl -i -nlp -e '/^(?:>)/ and s/.+/>$seq_name/' $id/$seq.fasta";
    }
    else {
        print " " x 4 . "Keep original names.\n";
    }

    system "perl $FindBin::Bin/../util/bp_genbank2gff3.pl $id/$seq.gb";

    move( "$seq.gb.gff", "$id/$seq.gff" );
    if ($pure_gff) {
        system "perl -i -nlp -e '/^\#\#FASTA/ and last' $id/$seq.gff";
    }

    if ( $rename and defined $row->[3] ) {
        my $seq_name = $row->[3];
        print " " x 4 . "Rename .fasta and .gff\n";
        move( "$id/$seq.fasta", "$id/$seq_name.fasta" );
        move( "$id/$seq.gff",   "$id/$seq_name.gff" );
    }

    print "\n\n";
}

close $csv_fh;

__END__

cd ~/data/alignment/yeast_genome
perl ~/Scripts/withncbi/util/batch_get_seq.pl -r -p -f yeast_name_seq.csv 2>&1 | tee yeast_name_seq.log
