#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
use Path::Tiny;
use Text::CSV_XS;
use Bio::Taxon;
use Bio::DB::Taxonomy;
use DateTime::Format::Natural;
use List::MoreUtils qw(any all uniq);

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use MyUtil;

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

gr_strains.pl

=head1 SYNOPSIS

    perl gr_strains.pl -o prok_strains.csv

    perl gr_strains.pl --euk -o euk_strains.csv

=cut

# running options
my $gr_dir = path( $Config->{path}{gr} )->stringify;    # genome report
my $td_dir = path( $Config->{path}{td} )->stringify;    # taxdmp

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'euk'    => \my $euk, # eukaryotes instead of prokaryotes
    'output|o=s' => \( my $strain_file = "prok_strains.csv" ),
) or Getopt::Long::HelpMessage(1);

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
my DBI $dbh = DBI->connect("DBI:CSV:");

if ( !$euk ) {

    # Chromosomes_RefSeq ==> Chromosomes
    $dbh->{csv_tables}->{t0} = {
        eol            => "\n",
        sep_char       => "\t",
        file           => "$gr_dir/prokaryotes.txt",
        skip_first_row => 1,
        quote_char     => '',
        col_names      => [
            qw{ Organism_Name TaxID BioProject_Accession BioProject_ID Group SubGroup Size GC
                Replicons WGS Scaffolds Genes Proteins Release_Date Modify_Date Status Center
                BioSample_Accession Assembly_Accession Reference FTP_Path Pubmed_ID Strain
                }
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
            qw{ Organism_Name TaxID BioProject_Accession BioProject_ID Group SubGroup Size GC
                Assembly_Accession
                Replicons WGS Scaffolds Genes Proteins Release_Date Modify_Date Status Center
                BioSample_Accession
                }
        ],
    };
}

#----------------------------#
# select columns only needed
#----------------------------#
{
    $stopwatch->block_message("Write summary");

    my $query = q{
        SELECT
            t0.TaxID,
            t0.Organism_Name,
            t0.BioProject_Accession,
            t0.Group,
            t0.SubGroup,
            t0.Size,
            t0.GC,
            t0.Replicons,
            t0.WGS,
            t0.Scaffolds,
            t0.Release_Date,
            t0.Assembly_Accession,
            t0.Status
        FROM t0
        WHERE 1 = 1
    };

    my DBI $header_sth = $dbh->prepare($query);
    $header_sth->execute;
    $header_sth->finish;

    # prepare output csv file
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, ">", $strain_file;
    my @cols_name = map { s/^t[01]\.//; $_ } @{ $header_sth->{'NAME'} };
    $csv->print(
        $csv_fh,
        [   @cols_name,
            qw{ species species_id genus genus_id species_member
                genus_species_member genus_strain_member }
        ]
    );

    my @strs = (
        q{ AND t0.Replicons <> '-'
            AND t0.Replicons <> ''
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
        my DBI $join_sth = $dbh->prepare( $query . $str );
        $join_sth->execute;
        while ( my @row = $join_sth->fetchrow_array ) {
            for my $item (@row) {
                $item = undef if ( $item eq '-' );
            }

            # only keep chromosomes
            if ( defined $row[7] ) {
                my @accs;
                for my $s ( split /;\s*/, $row[7] ) {
                    next unless $s =~ /chromosome/;
                    my @parts = split /:/, $s;
                    $parts[1] =~ s/\/.+//;
                    push @accs, $parts[1];

                }
                $row[7] = join "|", @accs;
            }

            # find each strains' species and genus
            my $taxon_id = $row[0];
            my $name     = $row[1];

            # dedup, make taxon_id unique
            next if grep { $_ == $row[0] } @taxon_ids;

            my Bio::Taxon $node = $taxon_db->get_taxon( -taxonid => $taxon_id );
            if ( !$node ) {
                warn "Can't find taxon for $name\n";
                next;
            }

            my $species = MyUtil::find_ancestor( $node, 'species' );
            if ($species) {
                push @row, ( $species->scientific_name, $species->id );
            }
            else {
                warn "Can't find species for $name\n";
                next;
            }

            my $genus = MyUtil::find_ancestor( $node, 'genus' );
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
    system "wc -l $_" for "$gr_dir/prokaryotes.txt", "$gr_dir/eukaryotes.txt", $strain_file;
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

__END__
