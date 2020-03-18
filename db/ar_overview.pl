#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use DBI;
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
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db_name  = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'output|o=s'   => \my $outfile,
) or Getopt::Long::HelpMessage(1);

unless ($outfile) {
    $outfile = $db_name . "_overview.xlsx";
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
$stopwatch->start_message("Overviews for $db_name...");

my $dbh = DBI->connect( "dbi:mysql:$db_name:$server", $username, $password )
    or die "Cannot connect to MySQL database at $db_name:$server";
my $to_xlsx = AlignDB::ToXLSX->new(
    dbh     => $dbh,
    outfile => $outfile,
);

#----------------------------------------------------------#
# worksheet -- strains
#----------------------------------------------------------#
my $strains = sub {
    my $sheet_name = 'strains';
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(0);

    my $sql_query = q{
        SELECT  *
        FROM ar
        WHERE 1 = 1
    };

    {    # header
        my @names = $to_xlsx->sql2names($sql_query);
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # write contents
            # species' member, chr_number and genus_member
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- species
#----------------------------------------------------------#
my $species = sub {
    my $sheet_name = 'species';
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(0);

    {    # write header
        my @names = qw{ genus_id genus species_id species species_member
            genus_species_member genus_strain_member code };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
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
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- group
#----------------------------------------------------------#
my $group = sub {
    my $sheet_name = 'group';
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(1);

    {    # write header
        my @names = qw{ group_name genus_member species_member strain_member };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # write contents
        my $query_name = 'superkingdom';
        my $sql_query  = q{
            SELECT
                superkingdom group_name,
                COUNT(DISTINCT genus_id) genus_member,
                COUNT(DISTINCT species_id) species_member,
                COUNT(DISTINCT taxonomy_id) strain_member
            FROM
                ar
            GROUP BY group_name
        };
        $to_xlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # write contents
        my $query_name = 'Euk subgroup';
        my $sql_query  = q{
            SELECT
                (CONCAT(`group`, ', ', subgroup)) group_name,
                COUNT(DISTINCT genus_id) genus_member,
                COUNT(DISTINCT species_id) species_member,
                COUNT(DISTINCT taxonomy_id) strain_member
            FROM ar
            WHERE superkingdom = 'Eukaryota'
            GROUP BY group_name
        };
        $to_xlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    {    # write contents
        my $query_name = 'Prok group';
        my $sql_query  = q{
            SELECT
                `group` group_name,
                COUNT(DISTINCT genus_id) genus_member,
                COUNT(DISTINCT species_id) species_member,
                COUNT(DISTINCT taxonomy_id) strain_member
            FROM ar
            WHERE superkingdom != 'Eukaryota'
            GROUP BY group_name
        };
        $to_xlsx->write_sql(
            $sheet,
            {   query_name => $query_name,
                sql_query  => $sql_query,
            }
        );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- euk_group
#----------------------------------------------------------#
my $euk_group = sub {
    my $sheet_name = 'euk_group';
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(0);

    {    # write header
        my @names = qw{ group_name genus_member species_member strain_member
            genome chromosome scaffold contig };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
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
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }
    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheets -- subgroup_XXX
#----------------------------------------------------------#
my $subgroup_query = sub {
    my $subgroup   = shift;
    my $sheet_name = "subgroup_$subgroup";
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(0);

    {    # write header
        my @names = qw{ genus genus_id species_member strain_member
            genome chromosome scaffold contig };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # write contents
        my $sql_query = qq{
            SELECT
                t0.genus,
                t0.genus_id,
                t0.species_member,
                t0.strain_member,
                t1.strain_member,
                t2.strain_member,
                t3.strain_member,
                t4.strain_member
            FROM
                (SELECT
                    genus,
                    genus_id,
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
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }
    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheets -- genus_XXX
#----------------------------------------------------------#
my $genus_query = sub {
    my $genus      = shift;
    my $sheet_name = "genus_$genus";
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(0);

    {    # write header
        my @names = qw{ species strain_member
            genome chromosome scaffold contig };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
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
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }
    print "Sheet [$sheet_name] has been generated.\n";
};

{
    &$strains;
    &$species;
    &$group;
    &$euk_group;
    $subgroup_query->("Ascomycetes");
    $subgroup_query->("Basidiomycetes");
    $genus_query->("Saccharomyces");
    $genus_query->("Dictyostelium");
    $genus_query->("Oryza");
}

$stopwatch->end_message;
exit;

__END__
