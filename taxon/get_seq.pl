#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Path::Tiny qw();
use LWP::Simple;

my $id  = shift // die "Provide a valid GenBank accession\n";
my $dir = shift || ".";

if ( !-d $dir ) {
    Path::Tiny::path($dir)->mkpath;
}

# Some .gb files only contain assembling information, so we need download both types.
# When network connection is awful, downloads may be incomplete.
my %suffix_of = ( gb => "gb", fasta => "fa" );

for my $type (keys %suffix_of) {
    my $url_template
        = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        . "db=nucleotide"
        . "&id=%s&rettype=%s"
        . "&retmode=text";

    my $url = sprintf $url_template, $id, $type;
    my $file = Path::Tiny::path( $dir, "$id.$suffix_of{$type}" )->absolute->stringify;

    my $rc = LWP::Simple::mirror( $url, $file );
    printf "* URL: %s\n" . "* LOCAL: %s\n", $url, $file;
    printf "* RC: %s\n\n", $rc;
}

__END__

=head1 NAME

get_seq.pl - retrieve nucleotide sequences from NCBI via eutils

=head1 SYNOPSIS

    perl get_seq.pl <Accession> [Output directory]
    
    perl get_seq.pl NC_012920 .

=cut
