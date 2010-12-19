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

my $chr_name = 'chr21';
my $all_set    = AlignDB::IntSpan->new("18635000-24650000,34051474-41801474");
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

my @rhs_types;
my $rhs_type_sth = $rhs_dbh->prepare(
    q{
    SELECT distinct hotspot_type
    FROM hotspot
    }
);
$rhs_type_sth->execute;
while ( my @row = $rhs_type_sth->fetchrow_array ) {
    my ($type) = @row;
    push @rhs_types, $type;
}

my %sets;
my $rhs_query_sth = $rhs_dbh->prepare(
    q{
    SELECT hotspot_runlist
    FROM hotspot
    where hotspot_type = ?
    }
);
for my $type (@rhs_types) {
    $sets{$type} = AlignDB::IntSpan->new;
    $rhs_query_sth->execute($type);
    while ( my @row = $rhs_query_sth->fetchrow_array ) {
        my ($runlist) = @row;
        $sets{$type}->merge($runlist);
    }

    $sets{$type} = $sets{$type}->pad($pad_length);
}
$rhs_query_sth->finish;

$sets{none_a} = $all_set->diff( $sets{Chimp}, $sets{AfAmer}, );
$sets{none_e} = $all_set->diff( $sets{Chimp}, $sets{EurAmer}, );
$sets{none_c} = $all_set->diff( $sets{Chimp}, $sets{ChnAmer}, );
$sets{none_s} = $all_set->diff( $sets{Chimp}, $sets{AfAmer}, $sets{EurAmer}, $sets{ChnAmer}, );
$sets{hotspot} = $all_set->diff($sets{none_s});

$sets{Chimp_oa} = AlignDB::IntSpan->new;
for my $set ( $sets{Chimp}->sets ) {
    if ( $sets{AfAmer}->intersect($set)->is_empty ) {
        $sets{Chimp_oa}->merge($set);
    }
}

$sets{Chimp_oe} = AlignDB::IntSpan->new;
for my $set ( $sets{Chimp}->sets ) {
    if ( $sets{EurAmer}->intersect($set)->is_empty ) {
        $sets{Chimp_oe}->merge($set);
    }
}

$sets{Chimp_oc} = AlignDB::IntSpan->new;
for my $set ( $sets{Chimp}->sets ) {
    if ( $sets{ChnAmer}->intersect($set)->is_empty ) {
        $sets{Chimp_oc}->merge($set);
    }
}

$sets{AfAmer_o} = AlignDB::IntSpan->new;
for my $set ( $sets{AfAmer}->sets ) {
    if ( $sets{Chimp}->intersect($set)->is_empty ) {
        $sets{AfAmer_o}->merge($set);
    }
}

$sets{EurAmer_o} = AlignDB::IntSpan->new;
for my $set ( $sets{EurAmer}->sets ) {
    if ( $sets{Chimp}->intersect($set)->is_empty ) {
        $sets{EurAmer_o}->merge($set);
    }
}

$sets{ChnAmer_o} = AlignDB::IntSpan->new;
for my $set ( $sets{ChnAmer}->sets ) {
    if ( $sets{Chimp}->intersect($set)->is_empty ) {
        $sets{ChnAmer_o}->merge($set);
    }
}


for my $key ( sort keys %sets ) {
    DumpFile("$chr_name.$key.yaml", $sets{$key}->runlist);
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
