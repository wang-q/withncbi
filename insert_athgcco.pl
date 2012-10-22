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
use AlignDB::Run;
use AlignDB::Window;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Multi;
use AlignDB::Ofg;

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

my $tag   = "gc";
my $style = "center";
my $noclean; # do not clean ofg tables

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{feature}{batch};

my $multi;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'p|port=i'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    't|tag=s'      => \$tag,
    'style=s'      => \$style,
    'parallel=i'   => \$parallel,
    'batch=i'      => \$batch_number,
    'multi'        => \$multi,
    'noclean'      => \$noclean,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

my $data_of = {
    gc => "$FindBin::Bin/athgcco/gene_conv.csv",
    co => "$FindBin::Bin/athgcco/crossover_10k.csv",
};
my $data_file = $data_of->{$tag};

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Update data of $db...");

my @jobs;
{
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    if ( !$noclean ) {
        print "Emptying tables...\n";

        # empty tables
        $obj->empty_table( 'ofg',   'with_window' );
        $obj->empty_table( 'ofgsw', 'with_window' );
    }

    my @align_ids = @{ $obj->get_align_ids };

    while ( scalar @align_ids ) {
        my @batching = splice @align_ids, 0, $batch_number;
        push @jobs, [@batching];
    }
}

#----------------------------------------------------------#
# read data
#----------------------------------------------------------#
my @all_data;
my %chr_data_set;
{
    open my $data_fh, '<', $data_file;
    <$data_fh>;    # omit head line
    while ( my $string = <$data_fh> ) {
        next unless defined $string;
        chomp $string;
        my ( $chr, $start, $end, $type )
            = ( split /\,/, $string )[ 0, 1, 2, 3 ];
        next unless $chr   =~ /^\d+$/;
        next unless $start =~ /^\d+$/;
        next unless $end   =~ /^\d+$/;
        if ( $start > $end ) {
            ( $start, $end ) = ( $end, $start );
        }
        next if $end - $start < 10;
        my $set = AlignDB::IntSpan->new("$start-$end");
        push @all_data,
            { chr => $chr, set => $set, tag => $tag, type => $type, };
        if ( !exists $chr_data_set{$chr} ) {
            $chr_data_set{$chr} = AlignDB::IntSpan->new;
        }
        $chr_data_set{$chr}->merge($set);
    }
    close $data_fh;
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my $job       = shift;
    my @align_ids = @$job;

    my $obj;
    if ( !$multi ) {
        $obj = AlignDB->new(
            mysql  => "$db:$server",
            user   => $username,
            passwd => $password,
        );
    }
    else {
        $obj = AlignDB::Multi->new(
            mysql  => "$db:$server",
            user   => $username,
            passwd => $password,
        );
    }
    AlignDB::Ofg->meta->apply($obj);
    $obj->style($style);

    $obj->insert_ofg( \@align_ids, \@all_data, \%chr_data_set );
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}
exit;

__END__

=head1 NAME

    insert_spo11.pl - Add annotation info to alignDB

=head1 SYNOPSIS

    insert_spo11.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --server          MySQL server IP/Domain name
       --db              database name
       --username        username
       --password        password

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
