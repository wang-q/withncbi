#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use Text::CSV_XS;

my $yml = LoadFile("DGRP.yml");

my $csv = Text::CSV_XS->new( { binary => 1 } )
    or die "Cannot use CSV: " . Text::CSV_XS->error_diag();
$csv->eol("\n");
open my $out_fh, ">", "dgrp.csv";
open my $ftp_fh, ">", "dgrp.ftp.txt";
for my $name ( sort keys %{$yml} ) {
    print "$name\n";

    for my $srx ( sort keys %{ $yml->{$name} } ) {
        my $info = $yml->{$name}{$srx};
        print " " x 4, "$srx\n";

        my $platform = $info->{platform};
        print " " x 8, $info->{platform}, "\n";

        my $layout = $info->{layout};
        print " " x 8, $layout, "\n";

        next unless $platform =~ /illumina|solexa/i;
        next unless $layout   =~ /pair/i;

        for my $i ( 0 .. scalar @{ $info->{srr} } - 1 ) {
            my $srr = $info->{srr}[$i];
            my $url = $info->{downloads}[$i];

            my $rg_str
                = '@RG'
                . "\tID:$srr"
                . "\tLB:$srx"
                . "\tPL:ILLUMINA"
                . "\tSM:$name";
            $csv->print( $out_fh,
                [ $name, $srx, $platform, $layout, $srr, $rg_str ] );
            print {$ftp_fh} $url, "\n";
        }
    }
}
close $ftp_fh;
close $out_fh;
