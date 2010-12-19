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

my $rhs_db = $Config->{rhs}->{db};

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
    'rhs=s'      => \$rhs_db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Update rhs-related tables of $db...");

my $obj = AlignDB::Ofg->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

# Database handler
my $dbh = $obj->dbh;

# rhs db handler
my $rhs_dbh
    = DBI->connect( "dbi:mysql:$rhs_db:$server", $username, $password );

# nature05 data is in chr21
my $chr_nature05 = 'chr21';

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

# empty tables
$obj->empty_ofg_tables;

# alignments
my $align_sth = $dbh->prepare(
    qq{
    SELECT t.align_id
    FROM target t, sequence s, chromosome c
    where t.seq_id = s.seq_id
    AND s.chr_id = c.chr_id
    AND c.chr_name = '$chr_nature05'
    ORDER BY align_id
    }
);
$align_sth->execute;
my $align_ids = $align_sth->selectcol_arrayref;
$align_sth->finish;

# all rhs
my @all_rhs;
my %chr_rhs_set;
$chr_rhs_set{$chr_nature05} = AlignDB::IntSpan->new;
my $rhs_query_sth = $rhs_dbh->prepare(
    q{
    SELECT hotspot_runlist,
           hotspot_type
      FROM hotspot
    }
);
$rhs_query_sth->execute;
while ( my @row = $rhs_query_sth->fetchrow_array ) {
    my ( $runlist, $tag ) = @row;
    my $set = AlignDB::IntSpan->new($runlist);
    push @all_rhs, { chr => $chr_nature05, set => $set, tag => $tag };
    $chr_rhs_set{$chr_nature05}->merge($set);
}
$rhs_query_sth->finish;

#----------------------------#
# INSERT INTO ofg and ofgsw
#----------------------------#
$obj->insert_ofg( $align_ids, \@all_rhs, \%chr_rhs_set );

$stopwatch->end_message;
exit;

__END__

=head1 NAME

    insert_gene.pl - Add annotation info to alignDB

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
