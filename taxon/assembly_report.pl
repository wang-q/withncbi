#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Path::Tiny qw();
use Text::CSV_XS;

my $file = shift || die "Need a NCBI ASSEMBLY report";

my @columns = (
    "Organism name",
    "Taxid",
    "Assembly name",
    "Infraspecific name",
    "BioSample",
    "BioProject",
    "Submitter",
    "Date",
    "Assembly type",
    "Release type",
    "Assembly level",
    "Genome representation",
    "WGS project",
    "Assembly method",
    "Genome coverage",
    "Sequencing technology",
    "RefSeq category",
    "RefSeq assembly accession",
    "GenBank assembly accession",
);

#@type Text::CSV_XS
my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );

if ( $file eq "header" ) {
    $csv->print( *STDOUT, [ map { s/\s+/_/g; $_ } @columns ] );
    exit;
}
elsif ( !Path::Tiny::path($file)->is_file ) {
    die "The input file [$file] doesn't exist";
}

#----------------------------#
# collect
#----------------------------#
my $content = Path::Tiny::path($file)->slurp;

my @c;
for my $key (@columns) {
    if ( $content =~ /$key\:\s*([\w .-]+)/m ) {
        push @c, $1;
    }
    else {
        push @c, "";
    }
}
$csv->print( *STDOUT, \@c );

__END__
