#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;

use IO::All;
use Text::CSV_XS;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_file;

my $strain_name;
my $taxon_id;

my $genbank;
my $ucsc;
my $scaffold;
my $nuclear;
my $chromosome;

my $man  = 0;
my $help = 0;

GetOptions(
    'help'      => \$help,
    'man'       => \$man,
    'f|file=s'  => \$in_file,
    'n|name=s'  => \$strain_name,
    't|taxon=i' => \$taxon_id,
    'genbank'   => \$genbank,
    'ucsc'      => \$ucsc,
    'scaffold'  => \$scaffold,
    'nuclear'   => \$nuclear,
    'chromosome'   => \$chromosome,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

die
    "Provide a input file (like GCA_000146045.2.assembly.txt) or a remote url.\n"
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
$csv->print( *STDOUT,
    [ '#strain_name', 'accession', 'strain_taxon_id', 'seq_name' ] );

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

    # strain_name,accession,strain_taxon_id,seq_name # This line is needed
    # S288c,NC_001133,559292,I
    my $accession = $genbank ? $fields[4] : $fields[6];
    $accession =~ s/\.\d+//;

    my $seq_name = $ucsc ? $fields[9] : $fields[0];

    $csv->print( *STDOUT, [ $strain_name, $accession, $taxon_id, $seq_name ] );
}

__END__

=head1 NAME

assemble_csv.pl - convert NCBI assemble report to a .csv file for batch_get_seq.pl

=head1 SYNOPSIS

    perl assemble_csv.pl <-f assembly.txt> [options]
      Options:
        --help              brief help message
        --man               full documentation
        -f, --file          input assemble report file
                            GCA_000149445.2.assembly.txt
                            ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt
        -n, --name          name of the strain, optional
        -t, --taxon         taxonomy id of the strain, optional
        --genbank           genbank instead of refseq
        --ucsc              UCSC-style sequence names
        --scaffold          include scaffolds
        --nuclear           exclude non-nuclear (Mt, Pt, Plasmids)

=head1 EXAMPLE

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt

=cut
