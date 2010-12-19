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

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

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
my $db       = $Config->{rhs}->{db};              # recombination hotspot

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
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Init $db...");

#----------------------------#
# call mysql
#----------------------------#
my $cmd    = "mysql -h$server -P$port -u$username -p$password ";
my $drop   = "-e \"DROP DATABASE IF EXISTS $db;\"";
my $create = "-e \"CREATE DATABASE $db;\"";
my $init   = $init_sql || "$FindBin::Bin/../rhs.sql";

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
my $excel = Spreadsheet::ParseExcel::Workbook->Parse('natgen05_modified.xls');

my $sheet_name_50k = "results50kb";
my $sheet_name_70k = "results70kb";
my $regex          = qr{chr(\d+)\:(\d+)\-(\d+)};

my $hotspot_file  = 'hotspots.chr_pos.txt.new';
my $coldspot_file = 'coldspots.chr_pos.txt.new';

#----------------------------------------------------------#
# natgen05 data, Ptak SE,
# check her homepage http://email.eva.mpg.de/~ptak/research.html#recmap
#----------------------------------------------------------#
# read in window info
print "Insert to window...\n";
foreach my $sheet ( @{ $excel->{Worksheet} } ) {
    next unless $sheet->{Name} eq $sheet_name_50k;

    print "Worksheet $sheet->{Name}\n";

    my $insert = $dbh->prepare(
        'INSERT INTO window (
            window_id, window_chr, window_start, window_end, window_runlist
        )
        VALUES (
            NULL, ?, ?, ?, ? 
        )'
    );

    foreach my $row ( 2 .. $sheet->{MaxRow} ) {
        my $string = $sheet->{Cells}[$row][19]{Val};    # column T
        next unless defined $string;
        $string =~ $regex;
        my ( $chr, $start, $end ) = ( $1, $2, $3 );
        my $set    = AlignDB::IntSpan->new("$start-$end");
        $insert->execute( $chr, $start, $end, $set->runlist );
    }
    $insert->finish;
}

# read in hotspots
print "Insert to hotspot...\n";
foreach my $sheet ( @{ $excel->{Worksheet} } ) {
    next unless $sheet->{Name} eq $sheet_name_70k;

    print "Worksheet $sheet->{Name}\n";

    my $insert = $dbh->prepare(
        'INSERT INTO hotspot (
            hotspot_id, window_id, hotspot_chr, hotspot_start,
            hotspot_end, hotspot_runlist, hotspot_type
        )
        VALUES (
            NULL, ?, ?, ?, 
            ?, ?, ?
        )'
    );

    foreach my $row ( 2 .. $sheet->{MaxRow} ) {
        my $string = $sheet->{Cells}[$row][15]{Val};    # column P
        next unless defined $string;
        insert_hotspot( $dbh, $string, 'Chimp' );
    }

    foreach my $row ( 2 .. $sheet->{MaxRow} ) {
        my $string = $sheet->{Cells}[$row][30]{Val};    # column AE
        next unless defined $string;
        insert_hotspot( $dbh, $string, 'AfAmer' );
    }

    foreach my $row ( 2 .. $sheet->{MaxRow} ) {
        my $string = $sheet->{Cells}[$row][45]{Val};    # column AT
        next unless defined $string;
        insert_hotspot( $dbh, $string, 'ChnAmer' );
    }

    foreach my $row ( 2 .. $sheet->{MaxRow} ) {
        my $string = $sheet->{Cells}[$row][60]{Val};    # column BI
        next unless defined $string;
        insert_hotspot( $dbh, $string, 'EurAmer' );
    }
}

#
#my $window_sth = $dbh->prepare(
#    "SElect w.window_id, w.window_runlist, h.hotspot_runlist
#    from window w
#    inner join hotspot h on w.window_id = h.window_id
#    where h.hotspot_type = ?"
#);
#print "\n" x 2;
#$window_sth->execute( 'Chimp' );
#to_csv($window_sth);
#print "\n" x 2;
#$window_sth->execute( 'AfAmer' );
#to_csv($window_sth);

#----------------------------------------------------------#
# science05 data, Myers S, see reference 8
#----------------------------------------------------------#

print "Insert to spot...\n";

# read in hotspot info
open my $hotspot_fh, '<', $hotspot_file;
while (my $string = <$hotspot_fh>) {
    next unless defined $string;
    insert_spot( $dbh, $string, 'hotspot' );
}
close $hotspot_fh;

# read in coldspot info
open my $coldspot_fh, '<', $coldspot_file;
while (my $string = <$coldspot_fh>) {
    next unless defined $string;
    insert_spot( $dbh, $string, 'coldspot' );
}
close $coldspot_fh;

$stopwatch->end_message();
exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#

sub find_window {
    my $dbh = shift;
    my $set = shift;

    my $window_sth = $dbh->prepare(
        'SElect window_id
        from window
        where window_start <= ?
        and window_end >= ?'
    );
    $window_sth->execute( $set->min, $set->max );
    my $tbl_ary_ref = $window_sth->fetchall_arrayref;
    warn "Can't find unique window\n" if @$tbl_ary_ref != 1;

    return $tbl_ary_ref->[0][0];
}

sub insert_hotspot {
    my $dbh        = shift;
    my $pos_string = shift;
    my $type       = shift;

    my $regex          = qr{(chr\d+)\:(\d+)\-(\d+)};
    my $insert = $dbh->prepare(
        'INSERT INTO hotspot (
            hotspot_id, window_id, hotspot_chr, hotspot_start,
            hotspot_end, hotspot_runlist, hotspot_type
        )
        VALUES (
            NULL, ?, ?, ?, 
            ?, ?, ?
        )'
    );

    $pos_string =~ $regex;
    my ( $chr, $start, $end ) = ( $1, $2, $3 );
    my $set = AlignDB::IntSpan->new("$start-$end");
    my $window_id = find_window( $dbh, $set );
    $insert->execute( $window_id, $chr, $start, $end, $set->runlist,
        $type );

    $insert->finish;

}

sub insert_spot {
    my $dbh        = shift;
    my $pos_string = shift;
    my $type       = shift;
    
    # chr6_cox_hap1:4071192-4075201 should be omitted
    my $regex          = qr{(chr[\dXY]+)(?:.*?)\:(\d+)\-(\d+)}; 
    my $insert = $dbh->prepare(
        q{
        INSERT INTO spot (
            spot_id, spot_chr, spot_start,
            spot_end, spot_runlist, spot_type
        )
        VALUES (
            NULL, ?, ?,
            ?, ?, ?
        )
        }
    );

    $pos_string =~ $regex;
    my ( $chr, $start, $end ) = ( $1, $2, $3 );
    my $set = AlignDB::IntSpan->new("$start-$end");
    $insert->execute( $chr, $start, $end, $set->runlist,
        $type );

    $insert->finish;

}

sub to_csv {
    my $sth = shift;

    require Text::CSV;
    my $csv = Text::CSV->new;

    # header line
    my @columns = @{ $sth->{NAME} };
    $csv->combine(@columns);
    print $csv->string, "\n";

    # all others
    while ( my @row = $sth->fetchrow_array ) {
        $csv->combine(@row);
        print $csv->string, "\n";
    }

}

__END__

=head1 NAME

    bootstrap.pl - bootstrap the indel-induced mutation rate increase

=head1 SYNOPSIS

    apply_sql.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --file            csv file name
       --rep             times of bootstrap simulations

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

