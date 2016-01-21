#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Roman;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all uniq zip);

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Window;
use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
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

my @tags;
my @types;
my $style = "center_intact";
my $noclean;    # do not clean ofg tables

my @files;

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $multi;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    's|server=s'   => \$server,
    'P|port=i'     => \$port,
    'd|db=s'       => \$db,
    'u|username=s' => \$username,
    'p|password=s' => \$password,
    'f|file=s'     => \@files,
    'tag=s'        => \@tags,
    'type=s'       => \@types,
    'style=s'      => \$style,
    'parallel=i'   => \$parallel,
    'batch=i'      => \$batch_number,
    'noclean'      => \$noclean,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

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
        $obj->create_column( "ofgsw", "ofgsw_dG", "DOUBLE" );
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
# prepare
my @args = zip @files, @tags, @types;

my @all_data;
my %chr_data_set;

while (@args) {
    my $file = shift @args;
    my $tag  = shift @args;
    my $type = shift @args;
    
    print "File [$file]\n";

    open my $data_fh, '<', $file;
    while ( my $string = <$data_fh> ) {
        next unless defined $string;
        chomp $string;
        my ( $chr, $start, $end )
            = ( split /\t/, $string )[ 0, 1, 2 ];
        next unless $chr =~ /^\w+$/;
        $chr =~ s/chr0?//i;
        next unless $start =~ /^\d+$/;
        next unless $end =~ /^\d+$/;
        if ( $start > $end ) {
            ( $start, $end ) = ( $end, $start );
        }
        my $set = AlignDB::IntSpan->new("$start-$end");
        push @all_data,
            { chr => $chr, set => $set, tag => $tag, type => $type };
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

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    AlignDB::Ofg->meta->apply($obj);
    $obj->style($style);
    $obj->max_out_distance(20);
    $obj->max_in_distance(20);

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

    insert_bed.pl - Add annotation info to alignDB

=head1 SYNOPSIS

perl init/init_alignDB.pl -d S288Cvsself
perl init/gen_alignDB_genome.pl -d S288Cvsself -t "4932,S288C" --dir d:\data\alignment\self_alignment\S288C\  --parallel 4
perl ofg/insert_bed.pl -d S288Cvsself --tag hot --type hot -f d:\wq\Scripts\alignDB\ofg\spo11\spo11_hot.bed --batch 10 --parallel 1

=cut
