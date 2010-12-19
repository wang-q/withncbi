#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Spreadsheet::ParseExcel;
use Statistics::Descriptive;
use DBI;

use Roman;

use AlignDB::IntSpan;

use FindBin;
use lib "$FindBin::Bin/../lib";
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

my $db = $Config->{gce}->{db};    # gene conversion events

my $init_sql;

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
    'init_sql=s' => \$init_sql,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Init $db...");

#----------------------------#
# call mysql
#----------------------------#
my $cmd    = "mysql -h$server -P$port -u$username -p$password ";
my $drop   = "-e \"DROP DATABASE IF EXISTS $db;\"";
my $create = "-e \"CREATE DATABASE $db;\"";
my $init   = $init_sql || "$FindBin::Bin/../gce.sql";

print "#drop\n" . "$cmd $drop\n";
system("$cmd $drop");
print "#create\n" . "$cmd $create\n";
system("$cmd $create");
print "#init\n" . "$cmd $db < $init\n";
system("$cmd $db < $init");

# Database handler
my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password );

print "\|$dbh\|\n";

#----------------------------------------------------------#
# Init
#----------------------------------------------------------#

my $gce_file = 'event_intervals.txt';
my $hs_file  = 'hot_spots.txt';

#----------------------------------------------------------#
# nature08 data, Nature 454, 479-485
#----------------------------------------------------------#

print "Insert to gce...\n";

# read in gce info
open my $gce_fh, '<', $gce_file;
while ( my $string = <$gce_fh> ) {
    next unless defined $string;
    chomp $string;
    insert_gce( $dbh, $string );
}
close $gce_fh;

print "Insert to gce hotspot...\n";

# read in hotspot info
open my $hs_fh, '<', $hs_file;
while ( my $string = <$hs_fh> ) {
    next unless defined $string;
    chomp $string;
    insert_hs( $dbh, $string );
}
close $hs_fh;

$stopwatch->end_message;
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#

sub insert_gce {
    my $dbh        = shift;
    my $pos_string = shift;

    my @fields = split /\t/, $pos_string;

    return unless $fields[0] =~ /chr(\d+)/;    # omit first line
    my $chr = $1;
    $chr = Roman($chr);    # ensembl use Roman number for S288c

    # chr first last size tetrad spores type
    my $insert = $dbh->prepare(
        q{
        INSERT INTO gce (
            gce_id, gce_chr, gce_start, gce_end, gce_runlist,
            gce_tetrad, gce_spores, gce_type
        )
        VALUES (
            NULL, ?, ?, ?, ?,
            ?, ?, ?
        )
        }
    );

    if ( $fields[1] =~ /Inf/i or $fields[2] =~ /Inf/i ) {
        print "Find Inf events\n";
        print $pos_string, "\n";
        return;
    }

    $fields[1] = int $fields[1];
    $fields[2] = int $fields[2];

    my $set = AlignDB::IntSpan->new("$fields[1]-$fields[2]");
    $insert->execute(
        $chr,       $fields[1], $fields[2], $set->runlist,
        $fields[4], $fields[5], $fields[6]
    );
    $insert->finish;

    return;
}

sub insert_hs {
    my $dbh        = shift;
    my $pos_string = shift;

    my @fields = split /\t/, $pos_string;

    return unless $fields[0] =~ /chr(\d+)/;    # omit first line
    my $chr = $1;
    $chr = Roman($chr);    # ensembl use Roman number for S288c

    # chr first last size tetrad spores type
    my $insert = $dbh->prepare(
        q{
        INSERT INTO hotspot (
            hotspot_id, hotspot_chr, hotspot_start, hotspot_end,
            hotspot_runlist, hotspot_type
        )
        VALUES (
            NULL, ?, ?, ?,
            ?, ?
        )
        }
    );

    if ( $fields[1] =~ /Inf/i or $fields[2] =~ /Inf/i ) {
        print "Find Inf events\n";
        print $pos_string, "\n";
        return;
    }

    $fields[1] = int $fields[1];
    $fields[2] = int $fields[2];

    my $set = AlignDB::IntSpan->new("$fields[1]-$fields[2]");
    $insert->execute( $chr, $fields[1], $fields[2], $set->runlist,
        $fields[3] );
    $insert->finish;

    return;
}

__END__

=head1 NAME

    gene_gce.pl - generate a gce database

=head1 SYNOPSIS

    apply_sql.pl [options]
     Options:
       --help            brief help message
       --man             full documentation

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

