#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Stopwatch;
use AlignDB::ToXLSX;

use FindBin;

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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db_name  = $Config->{database}{db};

# running options
my $outfile = "gr_overview.xlsx";

my $man  = 0;
my $help = 0;

GetOptions(
    'help'         => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'd|db=s'       => \$db_name,
    'o|output=s'   => \$outfile,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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
        ( $sheet, $sheet_row )
            = $to_xlsx->write_header_sql( $sheet_name, \%option );
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
        ( $sheet, $sheet_row )
            = $to_xlsx->write_header_direct( $sheet_name, \%option );
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
        ( $sheet, $sheet_row )
            = $to_xlsx->write_header_direct( $sheet_name, \%option );
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
        ( $sheet, $sheet_row )
            = $to_xlsx->write_header_direct( $sheet_name, \%option );
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

=head1 NAME

gr_overview.pl - Overviews for NCBI GENOME_REPORTS

=head1 SYNOPSIS

    perl gr_overview.pl --db gr

=cut

