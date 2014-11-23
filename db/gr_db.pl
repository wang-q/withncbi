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
use File::Slurp;
use List::MoreUtils qw(any all uniq);

use FindBin;

use AlignDB::Stopwatch;

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

my $init_sql = "$FindBin::Bin/../init.sql";

# running options
my $strain_file = "prok_strains.csv";
my $gr_dir      = $Config->{path}{gr};

# append euk
my $append;

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
    'init_sql=s'   => \$init_sql,
    'file=s'       => \$strain_file,
    'append'       => \$append,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Init genome report DB...");

#----------------------------#
# call mysql
#----------------------------#
if ( !$append ) {
    $stopwatch->block_message("Create DB skeleton for $db_name");

    my $drh = DBI->install_driver("mysql");    # Driver handle object
    $drh->func( 'dropdb',   $db_name, $server, $username, $password, 'admin' );
    $drh->func( 'createdb', $db_name, $server, $username, $password, 'admin' );
    $drh->func( 'reload',   $db_name, $server, $username, $password, 'admin' );

    my $dbh
        = DBI->connect( "dbi:mysql:$db_name:$server", $username, $password );
    open my $infh, '<', $init_sql;
    my $content = do { local $/; <$infh> };
    close $infh;
    my @statements = grep {/\w/} split /;/, $content;
    for (@statements) {
        $dbh->do($_) or die $dbh->errstr;
    }
}

my $dbh = DBI->connect( "dbi:mysql:$db_name:$server", $username, $password );

#----------------------------#
# Filling table strain
#----------------------------#
{
    $stopwatch->block_message("Loading $strain_file");

    my $load_sth = $dbh->prepare(
        qq{
        LOAD DATA LOCAL INFILE '$strain_file'
        INTO TABLE gr
        FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        }
    );
    $load_sth->execute;
    $load_sth->finish;
}

#----------------------------#
# Add code from GENOME_REPORTS
#----------------------------#
if ( !$append ) {
    $stopwatch->block_message("Add code");

    # csv should have the same column number as sql
    # so I add column code here
    $dbh->do(q{ ALTER TABLE gr ADD COLUMN code text });

    my @references;
    for my $file (
        "$gr_dir/prok_reference_genomes.txt",
        "$gr_dir/prok_representative_genomes.txt"
        )
    {
        my @lines = read_file( $file, chomp => 1, );
        push @references, @lines;
    }

    my %code_of;
    for my $line ( uniq(@references) ) {
        my ( $name, $code ) = ( split /\t/, $line )[ 2, 7 ];
        if ( exists $code_of{$name} ) {
            $code_of{$name} .= "/$code";
        }
        else {
            $code_of{$name} = $code;
        }
    }

    my $update_sth = $dbh->prepare(
        qq{
        UPDATE  gr
        SET     code = ?
        WHERE   organism_name = ?
        }
    );
    for my $name ( sort keys %code_of ) {
        $update_sth->execute( $code_of{$name}, $name );
    }
    $update_sth->finish;
}

#----------------------------#
# Count members
#----------------------------#
{
    $stopwatch->block_message("Update species and genus memberships");

    # find species contains multiply strains
    my %species_member_of;
    my $species_sth = $dbh->prepare(
        qq{
        SELECT species, count(taxonomy_id)
        FROM   gr
        WHERE 1 = 1
        GROUP BY species
        }
    );
    $species_sth->execute;
    while ( my @row = $species_sth->fetchrow_array ) {
        my ( $species, $count ) = @row;
        $species_member_of{$species} = $count;
    }
    $species_sth->finish;

    # find genus contains multiply species
    my %genus_member_of;
    my $genus_sth = $dbh->prepare(
        qq{
        SELECT genus, count(distinct species), count(taxonomy_id)
        FROM   gr
        WHERE 1 = 1
        GROUP BY genus
        }
    );
    $genus_sth->execute;
    while ( my @row = $genus_sth->fetchrow_array ) {
        my ( $genus, $species_count, $strain_count ) = @row;
        $genus_member_of{$genus} = [ $species_count, $strain_count ];
    }
    $genus_sth->finish;

    my $update_sth = $dbh->prepare(
        qq{
        UPDATE  gr
        SET     species_member = ?,
                genus_species_member = ?,
                genus_strain_member = ?
        WHERE   taxonomy_id = ?
        }
    );
    my $id_sth = $dbh->prepare(
        qq{
        SELECT taxonomy_id, species, genus
        FROM   gr
        }
    );
    $id_sth->execute;
    while ( my @row = $id_sth->fetchrow_array ) {
        my ( $taxonomy_id, $species, $genus ) = @row;
        $update_sth->execute(
            $species_member_of{$species},  $genus_member_of{$genus}->[0],
            $genus_member_of{$genus}->[1], $taxonomy_id
        );
    }
    $id_sth->finish;
    $update_sth->finish;
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

__END__

=head1 NAME

gr_db.pl

=head1 SYNOPSIS

perl gr_db.pl --file prok_strains.csv

perl gr_db.pl --append --file euk_strains.csv

=cut
