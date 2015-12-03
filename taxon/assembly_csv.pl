#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use IO::All;
use Text::CSV_XS;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

assemble_csv.pl - convert NCBI assemble report to a .csv file for batch_get_seq.pl

=head1 SYNOPSIS

    perl assemble_csv.pl <-f assembly.txt> [options]
      Options:
        --help      -?          brief help message

        --file      -f  STR     input assemble report file
                                GCA_000149445.2.assembly.txt
                                ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt
        --name      -n  STR     name of the strain, optional
        --taxon     -t  INT     taxonomy id of the strain, optional
        --genbank               genbank instead of refseq
        --ucsc                  UCSC-style sequence names
        --scaffold              include scaffolds
        --nuclear               exclude non-nuclear (Mt, Pt, Plasmids)
        --chromosome            keep chromosomes only
        --length        INT     skip sequences short than this

=head1 EXAMPLE

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt

=cut

GetOptions(
    'help|?'     => sub { HelpMessage(0) },
    'file|f=s'   => \my $in_file,
    'name|n=s'   => \my $strain_name,
    'taxon|t=i'  => \my $taxon_id,
    'genbank'    => \my $genbank,
    'ucsc'       => \my $ucsc,
    'scaffold'   => \my $scaffold,
    'nuclear'    => \my $nuclear,
    'chromosome' => \my $chromosome,
    'length=i'   => \my $length,
) or HelpMessage(1);

die "Provide a input file (like GCA_000146045.2.assembly.txt) or a remote url.\n"
    unless defined $in_file;

# 0 Sequence-Name
# 1 Sequence-Role
# 2 Assigned-Molecule
# 3 Assigned-Molecule-Location/Type
# 4 GenBank-Accn
# 5 Relationship
# 6 RefSeq-Accn
# 7 Assembly-Unit
# 8 Sequence-Length
# 9 UCSC-style-name

# IO::All
my $handle = io($in_file);

# Text::CSV_XS
my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );

# Header line
$csv->print( *STDOUT, [ '#strain_name', 'accession', 'strain_taxon_id', 'seq_name' ] );

while ( defined( my $line = $handle->getline ) ) {
    chomp $line;

    # meta info lines
    if ( $line =~ /^#/ ) {
        if ( $line =~ /Organism name\:\s*([\w -]+)/ ) {
            if ( !$strain_name ) {
                $strain_name = $1;
                $strain_name =~ s/\W+/_/g;
            }
        }
        if ( $line =~ /Taxid\:\s*(\d+)/ ) {
            if ( !$taxon_id ) {
                $taxon_id = $1;
            }
        }
        next;
    }

    # field lines
    my @fields = split /\t/, $line;

    # omit scaffolds
    if ( !$scaffold and $fields[1] ne 'assembled-molecule' ) {
        next;
    }

    # omit non-nuclear
    if ( $nuclear and $fields[7] eq 'non-nuclear' ) {
        next;
    }

    # keep chromosomes only
    if ( $chromosome and $fields[3] ne 'Chromosome' ) {
        next;
    }

    # skip short sequences
    if ( $length and $fields[8] < $length ) {
        next;
    }

    # strain_name,accession,strain_taxon_id,seq_name # This line is needed
    # S288c,NC_001133,559292,I
    my $accession = $genbank ? $fields[4] : $fields[6];
    $accession =~ s/\.\d+//;

    my $seq_name = $ucsc ? $fields[9] : $fields[0];
    $seq_name =~ s/\W+/_/g;

    $csv->print( *STDOUT, [ $strain_name, $accession, $taxon_id, $seq_name ] );
}

__END__
