#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;
use MCE::Flow Sereal => 1;

use List::MoreUtils::PP;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Window;

use App::RL::Common;

use lib "$FindBin::RealBin/../lib";
use AlignDB;
use AlignDB::Position;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

insert_position.pl - Insert positions to alignDB

=head1 SYNOPSIS

    perl insert_position.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --password  -p  STR     password
        --file      -f  @STR    position files
        --tag           @STR    tags
        --type          @STR    types
        --style         STR     ofg style, default is [center_intact]
        --batch         INT     number of positions in one child process, default is [500]
        --noclean               do not clean ofg tables

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
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
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
    'batch=i'      => \( my $batch_number = 500 ),
    'noclean'      => \( my $noclean ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Update data of $db...");

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
}

#----------------------------------------------------------#
# read data
#----------------------------------------------------------#
# prepare
my @args = List::MoreUtils::PP::mesh @files, @tags, @types;

$stopwatch->block_message("position information");
print YAML::Syck::Dump \@args;

unless (@args) {
    die "No positions to be processed\n";
}

my @all_data;
while (@args) {
    my $file = shift @args;
    my $tag  = shift @args;
    my $type = shift @args;

    print "File [$file]\n";

    open my $data_fh, '<', $file;
    while ( my $string = <$data_fh> ) {
        next unless defined $string;
        chomp $string;

        my $info = App::RL::Common::decode_header($string);
        next unless defined $info->{chr_name};
        $info->{tag}  = $tag;
        $info->{type} = $type;

        push @all_data, $info;
    }
    close $data_fh;
}

#----------------------------#
# Insert positions
#----------------------------#
$stopwatch->block_message("Insert positions");

my $worker_insert = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $wid = MCE->wid;
    print "* Process task [$chunk_id] by worker #$wid\n";

    my @data = @{$chunk_ref};

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    my DBI $dbh = $obj->dbh;
    my $pos_finder = AlignDB::Position->new( dbh => $dbh );

    # insert into ofg
    my DBI $ofg_insert_sth = $dbh->prepare(
        q{
        INSERT INTO ofg (
            ofg_id, window_id, ofg_tag, ofg_type
        )
        VALUES (
            NULL, ?, ?, ?
        )
        }
    );

    my %info_of;
    for my $item (@data) {
        my ( $align_id, $dummy )
            = @{ $obj->find_align( $item->{chr_name}, $item->{chr_start}, $item->{chr_end} ) };
        if ( !defined $align_id ) {
            warn " " x 4, "Can't find align for this position\n";
            warn YAML::Syck::Dump $item;
            next;
        }
        elsif ( defined $dummy ) {
            warn " " x 4, "Overlapped alignment in this position\n";
            warn YAML::Syck::Dump $item;
        }

        my $target_info = $obj->get_target_info($align_id);

        # target runlist
        my $target_set = AlignDB::IntSpan->new( $target_info->{seq_runlist} );

        # insert internal indels, that are, indels in target_set
        # indels in query_set is equal to spans of target_set minus one
        my $internal_indel_flag = 1;

        my $item_start = $pos_finder->at_align( $align_id, $item->{chr_start} );
        my $item_end   = $pos_finder->at_align( $align_id, $item->{chr_end} );
        next if $item_start > $item_end;

        my $item_set = AlignDB::IntSpan->new;
        $item_set->add_pair( $item_start, $item_end );
        $item_set = $item_set->intersect($target_set);

        # window
        my ($cur_window_id) = $obj->insert_window( $align_id, $item_set, $internal_indel_flag );

        # insert to table
        $ofg_insert_sth->execute( $cur_window_id, $item->{tag}, $item->{type} );

        $info_of{ $obj->last_insert_id } = {
            align_id => $align_id,
            set      => $item_set,
        };
    }

    printf " " x 4 . "Insert %d positions\n", scalar keys %info_of;
    MCE->gather(%info_of);
};

MCE::Flow::init {
    chunk_size  => $batch_number,
    max_workers => $parallel,
};
my %ofg_info_of = mce_flow $worker_insert, \@all_data;
MCE::Flow::finish;

#----------------------------#
# Insert sw
#----------------------------#
$stopwatch->block_message("Insert sliding windows");

my $worker_sw = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $wid = MCE->wid;
    print "* Process task [$chunk_id] by worker #$wid\n";

    my @ofg_ids = @{$chunk_ref};

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    my DBI $dbh = $obj->dbh;
    my $window_maker = AlignDB::Window->new(
        sw_size          => 100,
        max_out_distance => 20,
        max_in_distance  => 20,
    );

    # prepare ofgsw_insert
    my DBI $ofgsw_insert = $dbh->prepare(
        q{
        INSERT INTO ofgsw (
            ofgsw_id, window_id, ofg_id,
            ofgsw_type, ofgsw_distance
        )
        VALUES (
            NULL, ?, ?,
            ?, ?
        )
        }
    );

    for my $ofg_id (@ofg_ids) {
        my $align_id = $ofg_info_of{$ofg_id}->{align_id};
        my AlignDB::IntSpan $ofg_set = $ofg_info_of{$ofg_id}->{set};

        my $target_info = $obj->get_target_info($align_id);

        # target runlist
        my $target_set = AlignDB::IntSpan->new( $target_info->{seq_runlist} );

        # insert internal indels, that are, indels in target_set
        # indels in query_set is equal to spans of target_set minus one
        my $internal_indel_flag = 1;

        if ( $style =~ /^edge/ ) {

            # outside rsw
            my @out_rsw
                = $window_maker->outside_window( $target_set, $ofg_set->min, $ofg_set->max );

            for my $outside_rsw (@out_rsw) {
                my ($cur_window_id)
                    = $obj->insert_window( $align_id, $outside_rsw->{set}, $internal_indel_flag );

                $ofgsw_insert->execute( $cur_window_id, $ofg_id,
                    $outside_rsw->{type}, $outside_rsw->{distance},
                );
            }

            # inside rsw
            my @in_rsw = $window_maker->inside_window( $target_set, $ofg_set->min, $ofg_set->max );

            for my $inside_rsw (@in_rsw) {
                my ($cur_window_id)
                    = $obj->insert_window( $align_id, $inside_rsw->{set}, $internal_indel_flag );

                $ofgsw_insert->execute( $cur_window_id, $ofg_id,
                    $inside_rsw->{type}, $inside_rsw->{distance},
                );
            }

            if ( $style ne 'edge_only' ) {

                # inside rsw 2
                # rsw2 start from -90, so there will be no conflicts with rsw
                my @in_rsw2
                    = $window_maker->inside_window2( $target_set, $ofg_set->min, $ofg_set->max );

                for my $inside_rsw (@in_rsw2) {
                    my ($cur_window_id)
                        = $obj->insert_window( $align_id, $inside_rsw->{set},
                        $internal_indel_flag );

                    $ofgsw_insert->execute( $cur_window_id, $ofg_id,
                        $inside_rsw->{type}, $inside_rsw->{distance},
                    );
                }
            }
        }
        elsif ( $style eq 'center' ) {
            my @center_rsw
                = $window_maker->center_window( $target_set, $ofg_set->min, $ofg_set->max );

            for my $rsw (@center_rsw) {
                my ($cur_window_id)
                    = $obj->insert_window( $align_id, $rsw->{set}, $internal_indel_flag );

                $ofgsw_insert->execute( $cur_window_id, $ofg_id, $rsw->{type}, $rsw->{distance}, );
            }
        }
        elsif ( $style eq 'center_intact' ) {
            my @center_rsw
                = $window_maker->center_intact_window( $target_set, $ofg_set->min, $ofg_set->max );

            for my $rsw (@center_rsw) {
                my ($cur_window_id)
                    = $obj->insert_window( $align_id, $rsw->{set}, $internal_indel_flag );

                $ofgsw_insert->execute( $cur_window_id, $ofg_id, $rsw->{type}, $rsw->{distance}, );
            }
        }
    }
};

MCE::Flow::init {
    chunk_size  => $batch_number,
    max_workers => $parallel,
};
mce_flow $worker_sw, [ keys %ofg_info_of ];
MCE::Flow::finish;

my $report = sprintf "Ofg in files: [%s]\tOfg inserted: [%d]", scalar @all_data,
    scalar keys %ofg_info_of;
$stopwatch->block_message($report);

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
