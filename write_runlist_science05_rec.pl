#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::IntSpan;
use AlignDB::Position;
use AlignDB::Stopwatch;
use AlignDB::Window;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# Database init values
my $server   = $Config->{database}->{server};
my $port     = $Config->{database}->{port};
my $username = $Config->{database}->{username};
my $password = $Config->{database}->{password};
my $rhs_db   = $Config->{rhs}->{db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'username=s' => \$username,
    'password=s' => \$password,
    'rhs=s'      => \$rhs_db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

my $all_set    = AlignDB::IntSpan->new("1-242951149");
my $pad_length = 10_000;
my $outfile = 'runlist.txt';

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Write rhs-related tables...");

# rhs db handler
my $rhs_dbh
    = DBI->connect( "dbi:mysql:$rhs_db:$server", $username, $password );


#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

my %sets;
my $rhs_query_sth = $rhs_dbh->prepare(
    q{
    SELECT spot_runlist
    FROM spot
    where spot_type = ?
    and spot_chr = 'chr2'
    }
);
for my $type (qw{hotspot coldspot}) {
    $sets{$type} = AlignDB::IntSpan->new;
    $rhs_query_sth->execute($type);
    while ( my @row = $rhs_query_sth->fetchrow_array ) {
        my ($runlist) = @row;
        $sets{$type}->merge($runlist);
    }

    $sets{$type} = $sets{$type}->pad($pad_length);
}
$rhs_query_sth->finish;

$sets{normal} = $all_set->diff( values %sets );

for my $key ( sort keys %sets ) {
    DumpFile("$key.yaml", $sets{$key}->runlist);
}

$stopwatch->end_message();
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
