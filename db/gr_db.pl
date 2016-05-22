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
use List::MoreUtils::PP;

use AlignDB::Stopwatch;

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

gr_db.pl

=head1 SYNOPSIS

    # linux, mac
    perl gr_db.pl --db gr --file prok_strains.csv
    
    perl gr_db.pl --db gr --append --file euk_strains.csv
    
    # windows
    perl gr_db.pl --db gr --file prok_strains.csv --gr d:/data/NCBI/genomes/GENOME_REPORTS
    
    perl gr_db.pl --db gr --append --file euk_strains.csv --gr d:/data/NCBI/genomes/GENOME_REPORTS

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server      = $Config->{database}{server} ),
    'port|P=i'     => \( my $port        = $Config->{database}{port} ),
    'db|d=s'       => \( my $db_name     = $Config->{database}{db} ),
    'username|u=s' => \( my $username    = $Config->{database}{username} ),
    'password|p=s' => \( my $password    = $Config->{database}{password} ),
    'init_sql=s'   => \( my $init_sql    = "$FindBin::RealBin/../init.sql" ),
    'file=s'       => \( my $strain_file = "prok_strains.csv" ),
    'gr=s' => \( my $gr_dir = path( $Config->{path}{gr} )->stringify ),
    'append' => \my $append,    # append euk
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Init genome report DB...");

my $dsn
    = "dbi:mysql:"
    . "database="
    . $db_name
    . ";host="
    . $server
    . ";port="
    . $port;

#----------------------------#
# call mysql
#----------------------------#
if ( !$append ) {
    $stopwatch->block_message("Create DB skeleton for $db_name");

    my $drh = DBI->install_driver("mysql");    # Driver handle object
    $drh->func( 'dropdb',   $db_name, $server, $username, $password, 'admin' );
    $drh->func( 'createdb', $db_name, $server, $username, $password, 'admin' );

    # don't need this and crash under win32
    #$drh->func( 'reload',   $db_name, $server, $username, $password, 'admin' );

    my DBI $dbh = DBI->connect( $dsn, $username, $password );
    my $content = path($init_sql)->slurp;
    my @statements = grep {/\w/} split /;/, $content;
    for (@statements) {
        $dbh->do($_) or die $dbh->errstr;
    }
}

my DBI $dbh = DBI->connect( $dsn, $username, $password );

#----------------------------#
# Filling table strain
#----------------------------#
{
    $stopwatch->block_message("Loading $strain_file");

    my DBI $load_sth = $dbh->prepare(
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
        my @lines = path($file)->lines( { chomp => 1, } );
        push @references, @lines;
    }

    # MOD = Model organism
    # COM = Community selected reference
    my %code_of;
    for my $line ( List::MoreUtils::PP::uniq(@references) ) {
        my ( $name, $code ) = ( split /\t/, $line )[ 2, 7 ];

        my @codes = split /\//, $code;
        if ( exists $code_of{$name} ) {
            push @codes, split( /\//, $code_of{$name} );
        }
        @codes = grep { !/^(CLI|PHY)$/ } List::MoreUtils::PP::uniq(@codes);
        if ( grep { $_ eq 'MOD' } @codes ) {
            $code_of{$name} = 'MOD';
        }
        elsif ( grep { $_ eq 'COM' } @codes ) {
            $code_of{$name} = 'COM';
        }
        elsif ( grep { $_ eq 'FGS' } @codes ) {
            $code_of{$name} = 'FGS';
        }
        else {
            $code_of{$name} = join "/", @codes;
        }
    }

    my DBI $update_sth = $dbh->prepare(
        q{
        UPDATE  gr
        SET     code = ?
        WHERE   organism_name = ?
        }
    );
    for my $name ( sort keys %code_of ) {
        if ( length $code_of{$name} ) {
            $update_sth->execute( $code_of{$name}, $name );
        }
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
    my DBI $species_sth = $dbh->prepare(
        q{
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
    my DBI $genus_sth = $dbh->prepare(
        q{
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

    my DBI $update_sth = $dbh->prepare(
        q{
        UPDATE  gr
        SET     species_member = ?,
                genus_species_member = ?,
                genus_strain_member = ?
        WHERE   taxonomy_id = ?
        }
    );
    my DBI $id_sth = $dbh->prepare(
        q{
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
