#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $tdna_file = "$FindBin::Bin/tdna/test.txt";
my $output;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'f|file=s'   => \$tdna_file,
    'o|output=s' => \$output,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

=head1 NAME

    insert_tdna.pl - Add annotation info to alignDB

=head1 SYNOPSIS

    insert_tdna.pl [options]
      Options:
        --help            brief help message
        --man             full documentation
        --server          MySQL server IP/Domain name
        --db              database name
        --username        username
        --password        password
        --ensembl         ensembl database name

=cut
