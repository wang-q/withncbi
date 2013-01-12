#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long::Descriptive;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Text::CSV_XS;
use Text::Table;
use DBIx::XHTML_Table;

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");
my $conf_db = $Config->{database};

my ( $opt, $usage ) = describe_options(
    "usage: %c %o",
    [ 'help|h', 'display this message' ],
    [],
    ['Database init values'],
    [ 'server|s=s',   'MySQL IP/Domain', { default => $conf_db->{server} } ],
    [ 'port|P=i',     'MySQL port',      { default => $conf_db->{port} } ],
    [ 'username|u=s', 'username',        { default => $conf_db->{username} } ],
    [ 'password|p=s', 'password',        { default => $conf_db->{password} } ],
    [ 'db|d=s',       'database name',   { default => $conf_db->{db} } ],
    [],
    [ 'query|q=s',  'SQL statement', { default => "SELECT * FROM meta" } ],
    [ 'output|o=s', 'output filename' ],
    [   'type|t=s',
        'output style (csv, neat, table, box and html)',
        { default => "csv" }
    ],
);

$usage->die(
    {   pre_text =>
            "Write sql query results to a file, supporting multiple styles\n"
    }
) if $opt->{help};

$opt->{output} = "$opt->{db}.$opt->{type}" unless $opt->{output};

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Do stat for $opt->{db}...");

my $obj = AlignDB->new(
    mysql  => "$opt->{db}:$opt->{server}",
    user   => $opt->{username},
    passwd => $opt->{password},
);

# Database handler
my $dbh = $obj->dbh;

# Open output file
open my $out_fh, '>', $opt->{output}
    or die("Cannot open output file $opt->{output}");

# Execute sql query
my $sth = $dbh->prepare( $opt->{query} )
    or die $dbh->errstr;
$sth->execute
    or die $sth->errstr;

if ( $opt->{type} eq 'csv' ) {
    my $csv = Text::CSV_XS->new;

    # header line
    my @columns = @{ $sth->{NAME} };
    $csv->combine(@columns);
    print {$out_fh} $csv->string, "\n";

    # all others
    while ( my @row = $sth->fetchrow_array ) {
        $csv->combine(@row);
        print {$out_fh} $csv->string, "\n";
    }
}
elsif ( $opt->{type} eq 'neat' ) {

    # header line
    my @columns = @{ $sth->{NAME} };
    print {$out_fh} DBI::neat_list( \@columns ), "\n";

    # all others
    while ( my @row = $sth->fetchrow_array ) {
        print {$out_fh} DBI::neat_list( \@row ), "\n";
    }
}
elsif ( $opt->{type} eq 'box' or $opt->{type} eq 'table' ) {
    my $is_box = $opt->{type} eq 'box';

    # header line
    my @columns = map +{ title => $_, align_title => 'center' },
        @{ $sth->{NAME} };
    my $c = 0;
    splice @columns, $_ + $c++, 0, \' | ' for 1 .. $#columns;
    my @header_border = ( $is_box ? \' |' : () );
    my $table = Text::Table->new( @header_border, @columns, @header_border );

    # all others
    while ( my @row = $sth->fetchrow_array ) {
        $table->load( \@row );
    }
    my $rule = $table->rule(qw/- +/);
    my @rows_border = ( $is_box ? $rule : () );
    print {$out_fh} join '', @rows_border, $table->title, $rule, $table->body,
        @rows_border;
}
elsif ( $opt->{type} eq 'html' ) {
    my $columns = $sth->{NAME};
    my $rows    = $sth->fetchall_arrayref;

    my $table = DBIx::XHTML_Table->new( $rows, $columns );
    $table->modify( table => { border => 1, } );
    $table->modify(
        th => {
            style => {
                color      => '#a9b9a9',
                background => '#444444',
            }
        }
    );
    $table->modify(
        tr => { style => { background => [ '#bacaba', '#cbdbcb' ] } } );

    print {$out_fh} $table->output;
}
else {
    die "Unknown output style type!\n";
}

close $out_fh;

$stopwatch->end_message;
exit;

__END__

