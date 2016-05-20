#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Path::Tiny;
use LWP::Simple;

my $id  = shift || "NC_001284";
my $dir = shift || ".";

if ( !-d $dir ) {
    path($dir)->mkpath;
}

# Some .gb files only contain assembling infomation, so we need download both
# types.
# When network connection is aweful, downloads may be incomplete.
my @types = qw(gb fasta);

for my $type (@types) {
    my $url_template
        = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        . "db=nucleotide"
        . "&id=%s&rettype=%s"
        . "&retmode=text";

    my $url = sprintf $url_template, $id, $type;
    my $file = path( $dir, "$id.$type" )->absolute->stringify;

    my $rc = LWP::Simple::mirror( $url, $file );
    printf "* URL: %s\n" . "* LOCAL: %s\n", $url, $file;
    printf "* RC: %s\n\n", $rc;
}

__END__

=head1 NAME

get_seq.pl - retrieve nucleotide sequences from NCBI via eutils

=head1 SYNOPSIS

    perl get_seq.pl <Accession> [Output directory]
    
    perl get_seq.pl NC_001284 .

=cut
