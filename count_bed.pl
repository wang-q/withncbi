#!/usr/bin/perl
use strict;
use warnings;

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
use AlignDB::Position;

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

#my $style = "intersect";    # or superset, subset
my $noclean;    # do not clean previous counts

my @files;

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

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

    #'style=s'      => \$style,
    'parallel=i' => \$parallel,
    'batch=i'    => \$batch_number,
    'noclean'    => \$noclean,
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
        print "Create columns...\n";

        # empty tables
        $obj->create_column( "gsw", "gsw_bed_count", "INT" );
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
my $all_data     = {};
my $chr_data_set = {};

for my $file (@files) {
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
        if ( !exists $all_data->{$chr} ) {
            $all_data->{$chr} = [];
        }
        push @{ $all_data->{$chr} }, $set;

        if ( !exists $chr_data_set->{$chr} ) {
            $chr_data_set->{$chr} = AlignDB::IntSpan->new;
        }
        $chr_data_set->{$chr}->merge($set);
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

    my $dbh = $obj->dbh;
    my $pos_finder = AlignDB::Position->new( dbh => $dbh );

    # gsw
    my $fetch_gsw = $dbh->prepare(
        q{
        SELECT 
            g.gsw_id, w.window_start, w.window_end
        FROM
            gsw g
                INNER JOIN
            window w ON g.window_id = w.window_id
                AND w.align_id = ?
        }
    );

    # update sth
    my $bed_count_sth = $dbh->prepare(
        q{
        UPDATE gsw
        SET gsw_bed_count = ?
        WHERE gsw_id = ?
        }
    );

    # for each alignment
    for my $align_id (@align_ids) {
        my $target_info    = $obj->get_target_info($align_id);
        my $chr_name       = $target_info->{chr_name};
        my $chr_start      = $target_info->{chr_start};
        my $chr_end        = $target_info->{chr_end};
        my $align_length   = $target_info->{align_length};
        my $target_runlist = $target_info->{seq_runlist};

        next if $chr_name =~ /rand|un|contig|hap|scaf/i;

        $obj->process_message($align_id);

        $chr_name =~ s/chr0?//i;
        my $chr_set = AlignDB::IntSpan->new("$chr_start-$chr_end");

        # chr_ofg_set has intersect with chr_set
        #   ? there are ofgs in this alignmet
        #   : there is no ofg
        next unless exists $chr_data_set->{$chr_name};
        next if $chr_data_set->{$chr_name}->intersect($chr_set)->is_empty;

        # get sets within this alignment
        my @align_ofgs;
        for my $set ( @{ $all_data->{$chr_name} } ) {
            my $iset = $set->intersect($chr_set);
            if ( $iset->is_not_empty ) {
                print ' ' x 4, "Find ofg: $chr_name $iset\n";
                push @align_ofgs, $set;
            }
        }
        if ( scalar @align_ofgs == 0 ) {
            warn "Match wrong ofg positions\n";
            next;
        }

        # transform to align coordinate
        for my $set (@align_ofgs) {
            my $start = $pos_finder->at_align( $align_id, $set->min );
            my $end   = $pos_finder->at_align( $align_id, $set->max );
            next if $start > $end;
            $set = AlignDB::IntSpan->new("$start-$end");
        }

        $fetch_gsw->execute($align_id);
        while ( my ( $gsw_id, $window_start, $window_end, )
            = $fetch_gsw->fetchrow_array )
        {
            my $window_set = AlignDB::IntSpan->new("$window_start-$window_end");
            my @window_ofgs
                = grep { $window_set->intersect($_)->is_not_empty } @align_ofgs;
            $bed_count_sth->execute( scalar @window_ofgs, $gsw_id );
        }
    }
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

{    # set null to 0
    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );

    my $dbh           = $obj->dbh;
    my $bed_count_sth = $dbh->prepare(
        q{
        UPDATE gsw
        SET gsw_bed_count = 0
        WHERE gsw_bed_count IS NULL
        }
    );
    $bed_count_sth->execute;
}

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

    count_bed.pl - Add annotation info to alignDB

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
