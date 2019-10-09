#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Path::Tiny;
use Text::CSV_XS;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use MyUtil qw(find_ancestor find_group);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

ar_strains.pl

=head1 SYNOPSIS

    perl ar_strains.pl -o ar_strains.csv

    perl ar_strains.pl --genbank -o ar_strains_genbank.csv

=cut

# running options
my $ar_dir = path( $Config->{path}{ar} )->stringify;    # assembly report
my $td_dir = path( $Config->{path}{td} )->stringify;    # taxdmp

GetOptions(
    'help|?'  => sub { HelpMessage(0) },
    'genbank' => \my $genbank, # genbank instead of refseq
    'output|o=s' => \( my $strain_file = "ar_strains.csv" ),
) or HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Writing strains' summary...");

$stopwatch->block_message("Load ncbi taxdmp.");
my $taxon_db = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -directory => $td_dir,
    -nodesfile => "$td_dir/nodes.dmp",
    -namesfile => "$td_dir/names.dmp",
);

#----------------------------#
# load tab sep. txt files
#----------------------------#
$stopwatch->block_message("Load ncbi genome report and bioproject summary.");
my $dbh = DBI->connect("DBI:CSV:");

# Chromosomes_RefSeq ==> Chromosomes
$dbh->{csv_tables}->{t0} = {
    eol      => "\n",
    sep_char => "\t",
    file     => $genbank
    ? "$ar_dir/assembly_summary_genbank.txt"
    : "$ar_dir/assembly_summary_refseq.txt",
    skip_first_row => 1,
    quote_char     => '',
    col_names      => [
        qw{ assembly_accession bioproject biosample wgs_master
            refseq_category taxid species_taxid organism_name
            infraspecific_name isolate version_status assembly_level
            release_type genome_rep seq_rel_date asm_name submitter
            gbrs_paired_asm paired_asm_comp ftp_path }
    ],
};

#----------------------------#
# select columns only needed
#----------------------------#
{
    $stopwatch->block_message("Write summary");

    my $query = qq{
        SELECT 
            t0.TaxID,
            t0.organism_name,
            t0.bioproject,
            t0.assembly_accession,
            t0.wgs_master,
            t0.refseq_category,
            t0.assembly_level,
            t0.genome_rep,
            t0.seq_rel_date,
            t0.asm_name,
            t0.ftp_path
        FROM   t0
        WHERE 1 = 1
        AND t0.version_status = 'latest'
        AND t0.genome_rep = 'Full'
    };

    my $header_sth = $dbh->prepare($query);
    $header_sth->execute;
    $header_sth->finish;

    # prepare output csv file
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, ">", $strain_file or die "$strain_file: $!";
    my @cols_name = map { s/^t[01]\.//; $_ } @{ $header_sth->{'NAME'} };
    $csv->print(
        $csv_fh,
        [   @cols_name,
            qw{ superkingdom group subgroup species species_id genus genus_id species_member
                genus_species_member genus_strain_member }
        ]
    );

    my @strs = (
        q{ AND t0.assembly_level = 'Complete Genome'
            ORDER BY t0.seq_rel_date },
        q{ AND t0.assembly_level = 'Chromosome'
            ORDER BY t0.seq_rel_date },
        q{ AND t0.assembly_level = 'Scaffold'
            ORDER BY t0.seq_rel_date },
        q{ AND t0.assembly_level = 'Contig'
            ORDER BY t0.seq_rel_date },
    );
    my @taxon_ids;
    for my $str (@strs) {
        my $join_sth = $dbh->prepare( $query . $str );
        $join_sth->execute;
        while ( my @row = $join_sth->fetchrow_array ) {
            for my $item (@row) {
                $item = undef if ( $item eq '-' );
                $item = undef if ( $item eq 'na' );
            }

            # find each strains' species and genus
            my $taxon_id = $row[0];
            my $name     = $row[1];

            # dedup, make taxon_id unique
            next if grep { $_ == $row[0] } @taxon_ids;

            my $node = $taxon_db->get_taxon( -taxonid => $taxon_id );
            if ( !$node ) {
                warn "Can't find taxon for $name\n";
                next;
            }

            # superkingdom, group, subgroup
            push @row, ( find_group($node) );

            my $species = find_ancestor( $node, 'species' );
            if ($species) {
                push @row, ( $species->scientific_name, $species->id );
            }
            else {
                warn "Can't find species for $name\n";
                next;
            }

            my $genus = find_ancestor( $node, 'genus' );
            if ($genus) {
                push @row, ( $genus->scientific_name, $genus->id );
            }
            else {
                warn "Can't find genus for $name\n";
                next;
            }

            push @row, ( undef, undef, undef );    # member numbers

            # write a line
            push @taxon_ids, $row[0];
            $csv->print( $csv_fh, \@row );
        }
        $join_sth->finish;
    }
    close $csv_fh;
}

if ( $^O ne "Win32" ) {
    print "\n";
    system "wc -l $_"
        for "$ar_dir/assembly_summary_refseq.txt",
        "$ar_dir/assembly_summary_genbank.txt",
        $strain_file;
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

__END__
