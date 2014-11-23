#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Data::Table;
use Data::Table::Excel qw(tables2xlsx);

use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../config.ini");

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
    'help|?'       => \$help,
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

my $dbh = DBI->connect( "dbi:mysql:$db_name:$server", $username, $password );

#----------------------------#
# worksheet -- strains
#----------------------------#
my $t_strains = Data::Table::fromSQL(
    $dbh, q{
        SELECT * FROM gr WHERE 1 = 1
    }
);

#----------------------------#
# worksheet -- species
#----------------------------#
my $t_species = Data::Table::fromSQL(
    $dbh, q{
        SELECT genus_id, genus, species_id, species, AVG(genome_size),
                AVG(gc_content), species_member, genus_species_member,
                genus_strain_member, MAX(CHAR_LENGTH(code)) code
        FROM gr
        WHERE 1 = 1
        GROUP BY species
    }
);

#----------------------------#
# worksheet -- gc_checklist
#----------------------------#
my $t_gc_checklist = Data::Table::fromSQL(
    $dbh, q{
        SELECT genus_id, genus, species_id, species, AVG(genome_size),
                AVG(gc_content), COUNT(*) count, MAX(CHAR_LENGTH(code))
        FROM gr
        WHERE 1 = 1 AND species_member > 2
        GROUP BY species_id
        ORDER BY species
    }
);
push @{ $t_gc_checklist->{header} }, ( "table", "tree", "align", "xlsx" );

#----------------------------#
# worksheet -- gr_gc_checklist
#----------------------------#
my $t_gr_gc_checklist = Data::Table::fromSQL(
    $dbh, q{
        SELECT genus_id, genus, species_id, species, AVG(genome_size),
                AVG(gc_content), COUNT(*) count, MAX(CHAR_LENGTH(code))
        FROM gr
        WHERE 1 = 1 AND status = 'Complete' AND species_member > 2
        GROUP BY species_id HAVING count > 2
        ORDER BY species
    }
);
push @{ $t_gr_gc_checklist->{header} }, ( "table", "tree", "align", "xlsx" );

#----------------------------------------------------------#
# write file
#----------------------------------------------------------#
tables2xlsx(
    $outfile,
    [ $t_strains, $t_species, $t_gc_checklist, $t_gr_gc_checklist, ],
    [ "strains",  "species",  "gc_checklist",  "gr_gc_checklist", ]
);

$stopwatch->end_message;
exit;

__END__

=head1 NAME

gr_overview.pl

=head1 SYNOPSIS

perl gr_overview.pl

=cut
