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
my $db       = $Config->{database}{db};

my $tdna_file = "$FindBin::Bin/tdna/test.txt";

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
    'tdna=s'      => \$tdna_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Update tdna of $db...");

my $obj = AlignDB::Ofg->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
    max_in_distance => 0, # don't need inside windows
    max_out_distance => 10, # don't need inside windows
);

# Database handler
my $dbh = $obj->dbh;

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#

# empty tables
$obj->empty_ofg_tables;

# alignments
my $align_ids = $obj->get_align_ids;

# all gce
my @all_tdna;
my %chr_tdna_set;
{
    # read in gce info
    open my $tdna_fh, '<', $tdna_file;
    while ( my $string = <$tdna_fh> ) {
        next unless defined $string;
        chomp $string;
        my $pos_str = (split /\t/, $string)[1];
        next unless $pos_str;
        my ($chr, $pos) = split /:/, $pos_str;
        $chr =~ s/chr0?//i;
        $pos =~ s/^0+//;
        next unless $chr =~ /^\d+$/;
        my $set = AlignDB::IntSpan->new($pos);
        push @all_tdna, {chr => $chr, set => $set, tag => 'all'};
        if ( !exists $chr_tdna_set{$chr} ) {
            $chr_tdna_set{$chr} = AlignDB::IntSpan->new;
        }
        $chr_tdna_set{$chr}->merge($set);
    }
    close $tdna_fh;
}

DumpFile('tdna.yml', {tdna => \@all_tdna, chr => \%chr_tdna_set});

#----------------------------#
# INSERT INTO ofg and ofgsw
#----------------------------#
$obj->insert_ofg( $align_ids, \@all_tdna, \%chr_tdna_set );

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB::Ofg->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}
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
