#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

tdna2bed.pl - Convert tdna files from SALK to beds

=head1 SYNOPSIS

    perl tdna2bed.pl -i <tdna> [options]
      Options:
        --help      -?          brief help message
        --input     -i  STR     input file
        --output    -p  STR     output file, default is [$input.bed]

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'input|i=s'  => \my $tdna_file,
    'output|o=s' => \my $output,
) or HelpMessage(1);

unless ($output) {
    $output = "$tdna_file.bed";
}

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
open my $tdna_fh, '<', $tdna_file;
open my $out_fh,  '>', $output;
while ( my $string = <$tdna_fh> ) {
    next unless defined $string;
    chomp $string;
    my $pos_str = ( split /\t/, $string )[1];
    next unless $pos_str;
    my ( $chr, $pos ) = split /:/, $pos_str;
    $chr =~ s/chr0?//i;
    $pos =~ s/^0+//;
    next unless $chr =~ /^\d+$/;

    print {$out_fh} "$chr\t$pos\t$pos\n";
}
close $tdna_fh;
close $out_fh;

exit;

__END__
