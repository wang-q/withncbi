#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Roman;
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

my $data_file = "$FindBin::Bin/spo11/S288C_Spo11_Multimap.txt";

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
    'file=s'      => \$data_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Update data of $db...");

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
my @all_data;
my %chr_data_set;
{
    # read in gce info
    open my $data_fh, '<', $data_file;
    <$data_fh>; # omit head line
    while ( my $string = <$data_fh> ) {
        next unless defined $string;
        chomp $string;
        my ($chr, $start, $end) = (split /\t/, $string)[1, 3, 4];
        $chr =~ s/chr0?//i;
        next unless $chr =~ /^\d+$/;
        $chr = Roman($chr);
        next unless $chr =~ /^\w+$/;
        next unless $start =~ /^\d+$/;
        next unless $end =~ /^\d+$/;
        if ($start > $end) {
            ($start, $end) = ($end, $start);
        }
        my $set = AlignDB::IntSpan->new("$start-$end");
        push @all_data, {chr => $chr, set => $set, tag => 'all'};
        if ( !exists $chr_data_set{$chr} ) {
            $chr_data_set{$chr} = AlignDB::IntSpan->new;
        }
        $chr_data_set{$chr}->merge($set);
    }
    close $data_fh;
}

DumpFile('spo11.yml', {data => \@all_data, chr => \%chr_data_set});

#----------------------------#
# INSERT INTO ofg and ofgsw
#----------------------------#
$obj->insert_ofg( $align_ids, \@all_data, \%chr_data_set );

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
