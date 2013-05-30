#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Text::CSV_XS;

use FindBin;

my $file = shift;
die "Provide a csv file\n" unless defined $file and -e $file;

my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
open my $csv_fh, "<", $file or die "$file: $!";
$csv->getline($csv_fh);    # bypass title line

while ( my $row = $csv->getline($csv_fh) ) {
    my $id  = $row->[0];
    my $seq = $row->[1];

    print "id\t$id\tseq\t$seq\n";
    system "perl $FindBin::Bin/../util/get_seq.pl $seq $id";
    system "perl $FindBin::Bin/../util/genbank2gff.pl $id/$seq.gb > $id/$seq.gff";

}
close $csv_fh;
