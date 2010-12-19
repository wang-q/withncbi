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
use AlignDB::WriteExcel;
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
my $db       = $Config->{database}->{db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'username=s' => \$username,
    'password=s' => \$password,
    'db=s'       => \$db,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

my $pad_length = 5_000;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Write rhs-related tables...");

# alignDB object
my $alignDB_obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

my $write_obj = AlignDB::WriteExcel->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);

my $dbh = $alignDB_obj->dbh();

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

# find quartiles
my $quantiles;
{
    my $sql_query = q{
        SELECT gene_feature6
        FROM gene
        WHERE gene_feature6 IS NOT NULL
    };
    my %option = ( sql_query => $sql_query, );
    $quantiles = $write_obj->quantile_sql( \%option, 5 );
}
print Dump $quantiles;

my @rec_levels = (
    [ 1,   $quantiles->[0], $quantiles->[1] ],
    [ 2,   $quantiles->[1], $quantiles->[2] ],
    [ 3,   $quantiles->[2], $quantiles->[3] ],
    [ 4,   $quantiles->[3], $quantiles->[4] ],
    [ 5,   $quantiles->[4], $quantiles->[5] ],
    [ 'h', 1.5,           $quantiles->[5] ],
    [ 'c', $quantiles->[0], 0.82 ],
);

my @align_ids = @{ $alignDB_obj->get_align_ids };

# select all $goal chromosomes in this database
my @chrs = @{ $alignDB_obj->get_chrs('target') };
my %chr_set_of;
my $chr_level_of = {};
for my $chr (@chrs) {
    my ( $chr_id, $chr_name, $chr_length ) = @$chr;
    my $chr_set = AlignDB::IntSpan->new("1-$chr_length");
    $chr_set_of{$chr_name} = $chr_set;
    $chr_level_of->{$chr_name} = {};
    for (@rec_levels) {
        $chr_level_of->{$chr_name}{ $_->[0] } = AlignDB::IntSpan->new;
    }
}

for my $align_id (@align_ids) {
    my ( $chr_id, $chr_name, $chr_length )
        = $alignDB_obj->get_target_chr_info($align_id);

    my $sql_query = q{
        SELECT w.window_runlist
          FROM gene g, window w
         WHERE g.window_id = w.window_id
           AND g.gene_feature6 between ? and ?
           AND w.align_id = ?
    };

    my $sth = $dbh->prepare($sql_query);

    for my $level (@rec_levels) {
        my ( $order, $low_border, $high_border ) = @$level;
        $sth->execute( $low_border, $high_border, $align_id );

        while ( my @row = $sth->fetchrow_array ) {
            my ($runlist) = @row;
            $chr_level_of->{$chr_name}{$order}->merge($runlist);
        }
    }
}

for my $chr_name ( keys %{$chr_level_of} ) {

    for my $order ( keys %{ $chr_level_of->{$chr_name} } ) {
        my $set = $chr_level_of->{$chr_name}{$order};
        $set = $set->pad($pad_length);
        $set = $set->intersect( $chr_set_of{$chr_name} );
        DumpFile( "$order.$chr_name.yaml", $set->runlist );
    }
}

for my $chr_name ( keys %{$chr_level_of} ) {
    my $chr_set = $chr_set_of{$chr_name};
    my $set_h   = $chr_level_of->{$chr_name}{h};
    my $set_c   = $chr_level_of->{$chr_name}{c};
    my $set_n   = $chr_set->diff( $set_h, $set_c );
    DumpFile( "n.$chr_name.yaml", $set_n->runlist );

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
