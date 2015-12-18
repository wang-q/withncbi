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

gr_overview.pl - Overviews for NCBI GENOME_REPORTS

=head1 SYNOPSIS

    perl gr_overview.pl --db gr

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
            FROM gr
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
            FROM gr
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
        my @headers = qw{ genus_id genus species_id species
            avg_genome_size avg_gc species_member genus_species_member
            genus_strain_member code };
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
                    AVG(genome_size),
                    AVG(gc_content),
                    species_member,
                    genus_species_member,
                    genus_strain_member
            FROM    gr
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
# worksheet -- gc_checklist
#----------------------------------------------------------#
my $gc_checklist = sub {
    my $sheet_name = 'gc_checklist';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{
            genus_id genus species_id species avg_genome_size avg_gc count code
            table tree align xlsx
        };
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
# worksheet -- gr_gc_checklist
#----------------------------------------------------------#
my $gr_gc_checklist = sub {
    my $sheet_name = 'gr_gc_checklist';
    my $sheet;
    my ( $sheet_row, $sheet_col );

    {    # write header
        my @headers = qw{
            group genus_id species_id species avg_genome_size avg_gc count code
            table tree align xlsx
        };
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
            SELECT  `group` group_name,
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
            AND status like '%Complete%'
            AND species not like '%Candidatus%'
            GROUP BY species_id
            HAVING count > 2 AND species_code > 0
            ORDER BY group_name, species
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
    &$gc_checklist;
    &$gr_gc_checklist;
}

$stopwatch->end_message;
exit;

__END__
