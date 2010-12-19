#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use AlignDB::IntSpan;
use AlignDB::Window;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Ofg;
use AlignDB::Position;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $db       = $Config->{database}->{db};

my $gce_db = $Config->{gce}->{db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'gce=s'      => \$gce_db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Update gce-related tables of $db...");

my $obj = AlignDB::Ofg->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# gce db handler
my $gce_dbh
    = DBI->connect( "dbi:mysql:$gce_db:$server", $username, $password );

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

# empty tables
$obj->empty_ofg_tables;

# alignments
my @align_ids = @{ $obj->get_align_ids };

# all gce
my @all_gce;
my %chr_gce_set;
{
    my $gce_query_sth = $gce_dbh->prepare(
        q{
        SELECT gce_chr,
               gce_runlist,
               gce_type
          FROM gce
        }
    );
    $gce_query_sth->execute;
    while ( my @row = $gce_query_sth->fetchrow_array ) {
        my ( $chr, $runlist, $type ) = @row;
        my $set = AlignDB::IntSpan->new($runlist);
        push @all_gce,
            { chr => $chr, set => $set, tag => 'all', type => $type };
        if ( !exists $chr_gce_set{$chr} ) {
            $chr_gce_set{$chr} = AlignDB::IntSpan->new;
        }
        $chr_gce_set{$chr}->merge($set);
    }
    $gce_query_sth->finish;
}

# gce hotspot
my @all_hs;
my %chr_hs_set;
{
    my $hs_query_sth = $gce_dbh->prepare(
        q{
        SELECT hotspot_chr,
               hotspot_runlist,
               hotspot_type
          FROM hotspot
        }
    );
    $hs_query_sth->execute;
    while ( my @row = $hs_query_sth->fetchrow_array ) {
        my ( $chr, $runlist, $type ) = @row;
        my $set = AlignDB::IntSpan->new($runlist);
        push @all_hs,
            { chr => $chr, set => $set, tag => 'hotspot', type => $type };
        if ( !exists $chr_hs_set{$chr} ) {
            $chr_hs_set{$chr} = AlignDB::IntSpan->new;
        }
        $chr_hs_set{$chr}->merge($set);
    }
    $hs_query_sth->finish;
}

#----------------------------#
# INSERT INTO ofg and ofgsw
#----------------------------#
$obj->insert_ofg( \@align_ids, \@all_gce, \%chr_gce_set );
$obj->insert_ofg( \@align_ids, \@all_hs,  \%chr_hs_set );

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    insert_gce.pl - Add annotation info to alignDB

=head1 SYNOPSIS

    insert_gene.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password
       --ensembl         ensembl database name
       

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
