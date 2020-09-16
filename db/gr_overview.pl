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

gr_overview.pl - Overviews for NCBI GENOME_REPORTS

=head1 SYNOPSIS

    perl gr_overview.pl --db gr

=cut

GetOptions(
    'help|?'       => sub { Getopt::Long::HelpMessage(0) },
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
        FROM gr
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
        my @names = qw{
            genus_id genus species_id species avg_genome_size avg_gc species_member
            genus_species_member genus_strain_member code
        };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # write contents
         # species' member, chr_number and genus_member
        my $sql_query = q{
            SELECT  genus_id,
                    genus,
                    species_id,
                    species,
                    AVG(genome_size),
                    AVG(gc_content),
                    species_member,
                    genus_species_member,
                    genus_strain_member
            FROM    gr
            WHERE   1 = 1
            GROUP BY species
        };
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gc_checklist
#----------------------------------------------------------#
my $gc_checklist = sub {
    my $sheet_name = 'gc_checklist';
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(0);

    {    # write header
        my @names = qw{
            genus_id genus species_id species avg_genome_size avg_gc count code
            table tree align xlsx
        };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # write contents
        my $sql_query = q{
            SELECT  genus_id,
                    genus,
                    species_id,
                    species,
                    AVG(genome_size),
                    AVG(gc_content),
                    COUNT(*) count,
                    MAX(CHAR_LENGTH(code))
            FROM gr
            WHERE   1 = 1
            AND species_member > 2
            GROUP BY species_id
            ORDER BY species
        };
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

#----------------------------------------------------------#
# worksheet -- gr_gc_checklist
#----------------------------------------------------------#
my $gr_gc_checklist = sub {
    my $sheet_name = 'gr_gc_checklist';
    my $sheet;
    $to_xlsx->row(0);
    $to_xlsx->column(0);

    {    # write header
        my @names = qw{
            group genus_id species_id species avg_genome_size avg_gc count code
            table tree align xlsx
        };
        $sheet = $to_xlsx->write_header( $sheet_name, { header => \@names } );
    }

    {    # write contents
        my $sql_query = q{
            SELECT  `subgroup` subgroup,
                    genus_id,
                    species_id,
                    species,
                    AVG(genome_size),
                    AVG(gc_content),
                    COUNT(*) count,
                    MAX(CHAR_LENGTH(code)) species_code
            FROM gr
            WHERE   1 = 1
            AND species_member > 2
            AND status like "%Complete%"
            AND species not like "%Candidatus%"
            GROUP BY species_id
            HAVING count > 2 AND species_code > 0
            ORDER BY subgroup, species
        };
        $to_xlsx->write_sql( $sheet, { sql_query => $sql_query, } );
    }

    print "Sheet [$sheet_name] has been generated.\n";
};

{
    &$strains;
    &$species;
    &$gc_checklist;
    &$gr_gc_checklist;
}

$stopwatch->end_message;
exit;

__END__
