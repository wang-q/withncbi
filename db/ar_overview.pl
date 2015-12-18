#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Stopwatch;
use AlignDB::ToXLSX;

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

ar_overview_tx.pl - Overviews for NCBI ASSEMBLY_REPORTS

=head1 SYNOPSIS

    perl ar_overview.pl --db ar_refseq
    perl ar_overview.pl --db ar_genbank -o ar_genebank_overview.xlsx

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db_name  = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'output|o=s'   => \my $outfile,
) or HelpMessage(1);

unless ($outfile) {
    $outfile = $db_name . "_overview.xlsx";
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
$stopwatch->start_message("Overviews for $db_name...");

my $to_xlsx = AlignDB::ToXLSX->new(
    mysql   => "$db_name:$server",
    user    => $username,
    passwd  => $password,
    outfile => $outfile,
);

#----------------------------------------------------------#
# worksheet -- strains
#----------------------------------------------------------#
my $strains = sub {
    my $sheet_name = 'strains';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my $sql_query = q{
            SELECT  *
            FROM ar
            WHERE 1 = 1
        };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ( $sheet, $sheet_row ) = $to_xlsx->write_header_sql( $sheet_name, \%option );
    }

    {    # write contents
            # species' member, chr_number and genus_member
        my $sql_query = q{
            SELECT  *
            FROM ar
            WHERE 1 = 1
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- species
#----------------------------------------------------------#
my $species = sub {
    my $sheet_name = 'species';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{ genus_id genus species_id species species_member
            genus_species_member genus_strain_member code };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row ) = $to_xlsx->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
            # species' member, chr_number and genus_member
        my $sql_query = q{
            SELECT  genus_id,
                    genus,
                    species_id,
                    species,
                    species_member,
                    genus_species_member,
                    genus_strain_member
            FROM    ar
            WHERE   1 = 1
            GROUP BY species
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- group
#----------------------------------------------------------#
my $group = sub {
    my $sheet_name = 'group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{ group_name genus_member species_member strain_member };
        ( $sheet_row, $sheet_col ) = ( 0, 1 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row ) = $to_xlsx->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $query_name = 'superkingdom';
        $sheet_row++;
        my $sql_query = q{
            SELECT 
                superkingdom group_name,
                COUNT(DISTINCT genus_id) genus_member,
                COUNT(DISTINCT species_id) species_member,
                COUNT(DISTINCT taxonomy_id) strain_member
            FROM
                ar
            GROUP BY group_name
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'Euk subgroup';
        $sheet_row++;
        my $sql_query = q{
            SELECT 
                (CONCAT(`group`, ', ', subgroup)) group_name,
                COUNT(DISTINCT genus_id) genus_member,
                COUNT(DISTINCT species_id) species_member,
                COUNT(DISTINCT taxonomy_id) strain_member
            FROM ar
            WHERE superkingdom = 'Eukaryota'
            GROUP BY group_name
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }

    {    # write contents
        my $query_name = 'Prok group';
        $sheet_row++;
        my $sql_query = q{
            SELECT 
                `group` group_name,
                COUNT(DISTINCT genus_id) genus_member,
                COUNT(DISTINCT species_id) species_member,
                COUNT(DISTINCT taxonomy_id) strain_member
            FROM ar
            WHERE superkingdom != 'Eukaryota'
            GROUP BY group_name
        };
        my %option = (
            query_name => $query_name,
            sql_query  => $sql_query,
            sheet_row  => $sheet_row,
            sheet_col  => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }

    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- euk_group
#----------------------------------------------------------#
my $euk_group = sub {
    my $sheet_name = 'euk_group';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{ group_name genus_member species_member strain_member
            genome chromosome scaffold contig };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row ) = $to_xlsx->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $sql_query = q{
            SELECT 
                t0.group_name,
                t0.genus_member,
                t0.species_member,
                t0.strain_member,
                t1.strain_member,
                t2.strain_member,
                t3.strain_member,
                t4.strain_member
            FROM
                (SELECT 
                    (CONCAT(`group`, ', ', subgroup)) group_name,
                        COUNT(DISTINCT genus_id) genus_member,
                        COUNT(DISTINCT species_id) species_member,
                        COUNT(taxonomy_id) strain_member
                FROM ar
                WHERE superkingdom = 'Eukaryota'
                GROUP BY group_name) t0
                    LEFT JOIN
                (SELECT 
                    (CONCAT(`group`, ', ', subgroup)) group_name,
                        COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE superkingdom = 'Eukaryota'
                        AND assembly_level LIKE '%Genome%'
                GROUP BY group_name) t1 ON t0.group_name = t1.group_name
                    LEFT JOIN
                (SELECT 
                    (CONCAT(`group`, ', ', subgroup)) group_name,
                        COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE superkingdom = 'Eukaryota'
                        AND assembly_level LIKE '%Chromosome%'
                GROUP BY group_name) t2 ON t0.group_name = t2.group_name
                    LEFT JOIN
                (SELECT 
                    (CONCAT(`group`, ', ', subgroup)) group_name,
                        COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE superkingdom = 'Eukaryota'
                        AND assembly_level = 'Scaffold'
                GROUP BY group_name) t3 ON t0.group_name = t3.group_name
                    LEFT JOIN
                (SELECT 
                    (CONCAT(`group`, ', ', subgroup)) group_name,
                        COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE superkingdom = 'Eukaryota'
                        AND assembly_level = 'Contig'
                GROUP BY group_name) t4 ON t0.group_name = t4.group_name
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }
    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheets -- subgroup_XXX
#----------------------------------------------------------#
my $subgroup_query = sub {
    my $subgroup   = shift;
    my $sheet_name = "subgroup_$subgroup";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{ genus species_member strain_member
            genome chromosome scaffold contig };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row ) = $to_xlsx->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $sql_query = qq{
            SELECT 
                t0.genus,
                t0.species_member,
                t0.strain_member,
                t1.strain_member,
                t2.strain_member,
                t3.strain_member,
                t4.strain_member
            FROM
                (SELECT 
                    genus,
                    COUNT(DISTINCT species_id) species_member,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE subgroup = '$subgroup'
                GROUP BY genus) t0
                    LEFT JOIN
                (SELECT 
                    genus,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE subgroup = '$subgroup'
                        AND assembly_level LIKE '%Genome%'
                GROUP BY genus) t1 ON t0.genus = t1.genus
                    LEFT JOIN
                (SELECT 
                    genus,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE subgroup = '$subgroup'
                        AND assembly_level LIKE '%Chromosome%'
                GROUP BY genus) t2 ON t0.genus = t2.genus
                    LEFT JOIN
                (SELECT 
                    genus,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE subgroup = '$subgroup'
                        AND assembly_level = 'Scaffold'
                GROUP BY genus) t3 ON t0.genus = t3.genus
                    LEFT JOIN
                (SELECT 
                    genus,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE subgroup = '$subgroup'
                        AND assembly_level = 'Contig'
                GROUP BY genus) t4 ON t0.genus = t4.genus
            ORDER BY t0.strain_member DESC
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }
    print "Sheet \"$sheet_name\" has been generated.\n";
};

#----------------------------------------------------------#
# worksheets -- genus_XXX
#----------------------------------------------------------#
my $genus_query = sub {
    my $genus      = shift;
    my $sheet_name = "genus_$genus";
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{ species strain_member
            genome chromosome scaffold contig };
        ( $sheet_row, $sheet_col ) = ( 0, 0 );
        my %option = (
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
            header    => \@headers,
        );
        ( $sheet, $sheet_row ) = $to_xlsx->write_header_direct( $sheet_name, \%option );
    }

    {    # write contents
        my $sql_query = qq{
            SELECT 
                t0.species,
                t0.strain_member,
                t1.strain_member,
                t2.strain_member,
                t3.strain_member,
                t4.strain_member
            FROM
                (SELECT 
                    species,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE genus = '$genus'
                GROUP BY species) t0
                    LEFT JOIN
                (SELECT 
                    species,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE genus = '$genus'
                        AND assembly_level LIKE '%Genome%'
                GROUP BY species) t1 ON t0.species = t1.species
                    LEFT JOIN
                (SELECT 
                    species,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE genus = '$genus'
                        AND assembly_level LIKE '%Chromosome%'
                GROUP BY species) t2 ON t0.species = t2.species
                    LEFT JOIN
                (SELECT 
                    species,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE genus = '$genus'
                        AND assembly_level = 'Scaffold'
                GROUP BY species) t3 ON t0.species = t3.species
                    LEFT JOIN
                (SELECT 
                    species,
                    COUNT(DISTINCT taxonomy_id) strain_member
                FROM ar
                WHERE genus = '$genus'
                        AND assembly_level = 'Contig'
                GROUP BY species) t4 ON t0.species = t4.species
            ORDER BY t0.strain_member DESC
        };
        my %option = (
            sql_query => $sql_query,
            sheet_row => $sheet_row,
            sheet_col => $sheet_col,
        );
        ($sheet_row) = $to_xlsx->write_content_direct( $sheet, \%option );
    }
    print "Sheet \"$sheet_name\" has been generated.\n";
};

{
    &$strains;
    &$species;
    &$group;
    &$euk_group;
    $subgroup_query->("Ascomycetes");
    $genus_query->("Saccharomyces");
    $genus_query->("Dictyostelium");
    $genus_query->("Aspergillus");
    $genus_query->("Candida");
    $genus_query->("Oryza");
}

$stopwatch->end_message;
exit;

__END__
