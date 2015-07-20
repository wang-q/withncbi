#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Text::CSV_XS;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use MyUtil qw(replace_home find_ancestor);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

# running options
my $gr_dir = replace_home( $Config->{path}{gr} );    # genome report
my $td_dir = replace_home( $Config->{path}{td} );    # taxdmp

# eukaryotes instead of prokaryotes
my $euk;

my $strain_file = "prok_strains.csv";

my $man  = 0;
my $help = 0;

GetOptions(
    'help'       => \$help,
    'man'        => \$man,
    'euk'        => \$euk,
    'o|output=s' => \$strain_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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

if ( !$euk ) {

    # Chromosomes_RefSeq ==> Chromosomes
    $dbh->{csv_tables}->{t0} = {
        eol            => "\n",
        sep_char       => "\t",
        file           => "$gr_dir/prokaryotes.txt",
        skip_first_row => 1,
        quote_char     => '',
        col_names      => [
            qw{ Organism_Name TaxID BioProject_Accession BioProject_ID Group
                SubGroup Size GC Chromosomes Chromosomes_INSDC Plasmids_RefSeq
                Plasmids_INSDC WGS Scaffolds Genes Proteins Release_Date Modify_Date
                Status Center }
        ],
    };
}
else {
    $dbh->{csv_tables}->{t0} = {
        eol            => "\n",
        sep_char       => "\t",
        file           => "$gr_dir/eukaryotes.txt",
        skip_first_row => 1,
        quote_char     => '',
        col_names      => [
            qw{ Organism_Name TaxID BioProject_Accession BioProject_ID Group
                SubGroup Size GC Assembly Chromosomes Organelles Plasmids WGS
                Scaffolds Genes Proteins Release_Date Modify_Date Status Center }
        ],
    };
}

#----------------------------#
# select columns only needed
#----------------------------#
{
    $stopwatch->block_message("Write summary");

    my $query = qq{
        SELECT 
            t0.TaxID,
            t0.Organism_Name,
            t0.BioProject_Accession,
            t0.Group,
            t0.SubGroup,
            t0.Size,
            t0.GC,
            t0.Chromosomes,
            t0.WGS,
            t0.Scaffolds,
            t0.Release_Date,
            t0.Status
        FROM t0
        WHERE 1 = 1
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
            qw{ species species_id genus genus_id species_member
                genus_species_member genus_strain_member }
        ]
    );

    my @strs = (

        q{ AND t0.Chromosomes <> '-'
            AND t0.Chromosomes <> ''
            ORDER BY t0.Release_Date },
        q{ AND t0.Status = 'Scaffold'
            AND t0.WGS <> '-'
            AND t0.WGS <> ''
            ORDER BY t0.Release_Date },
        q{ AND t0.Status = 'Contig'
            AND t0.WGS <> '-'
            AND t0.WGS <> ''
            ORDER BY t0.Release_Date },
    );
    my @taxon_ids;
    for my $str (@strs) {
        my $join_sth = $dbh->prepare( $query . $str );
        $join_sth->execute;
        while ( my @row = $join_sth->fetchrow_array ) {
            for my $item (@row) {
                $item = undef if ( $item eq '-' );
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
        for "$gr_dir/prokaryotes.txt", "$gr_dir/eukaryotes.txt",
        $strain_file;
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

__END__

=head1 NAME

    gr_strains.pl

=head1 SYNOPSIS

    perl gr_strains.pl -o prok_strains.csv

    perl gr_strains.pl --euk -o euk_strains.csv

=cut

