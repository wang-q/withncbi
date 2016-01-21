#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Roman;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all uniq zip);

use AlignDB::IntSpan;
use AlignDB::Run;
use AlignDB::Window;
use AlignDB::Stopwatch;

use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Ofg;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

insert_bed.pl - Insert beds to alignDB

=head1 SYNOPSIS

    perl insert_bed.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --password  -p  STR     password
        --file      -f  @STR    bed files
        --tag           @STR    tags
        --type          @STR    types
        --style         STR     ofg style, default is [center_intact]
        --batch         INT     number of alignments in one child process
        --noclean               do not clean ofg tables
        --dG                    calculate deltaG

        
    Styles:
        edge
        ------+------------------------------------+--------
           2 1 -1 -2     -89 -90  -90 -89     -2 -1 1 2
        
        edge_only
        ------+-----------------------+--------
           2 1 -1 -2             -2 -1 1 2
        
        center
        ------+-----------------+--------
        7 6 5 4 3 2 1 0 1 2 3 4 5 6 7 8
        
        center_intact
        ------+-----------------+--------
         3 2 1        0          1 2 3 4

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server       = $Config->{database}{server} ),
    'port|P=i'     => \( my $port         = $Config->{database}{port} ),
    'db|d=s'       => \( my $db           = $Config->{database}{db} ),
    'username|u=s' => \( my $username     = $Config->{database}{username} ),
    'password|p=s' => \( my $password     = $Config->{database}{password} ),
    'file|f=s'     => \my @files,
    'tag=s'        => \my @tags,
    'type=s'       => \my @types,
    'style=s'      => \( my $style        = "center_intact" ),
    'parallel=i'   => \( my $parallel     = $Config->{generate}{parallel} ),
    'batch=i'      => \( my $batch_number = $Config->{generate}{batch} ),
    'noclean'      => \( my $noclean ),
    'dG'           => \( my $deltaG ),
) or HelpMessage(1);

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

$stopwatch->block_message("bed information");
print Dump \@args;

unless (@args) {
    die "No bed to be processed\n";    
}

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
        push @all_data, { chr => $chr, set => $set, tag => $tag, type => $type };
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

    if ($deltaG) {
        $obj->insert_dG(1);
    }
    $obj->insert_ofg( \@align_ids, \@all_data, \%chr_data_set );
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
$stopwatch->block_message("Start update");

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
