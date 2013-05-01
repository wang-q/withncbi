#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use File::Find::Rule;
use List::MoreUtils qw(any all uniq zip);
use Set::Scalar;

use Text::CSV_XS;

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

# executable file location
my $dir_result    = "/home/wangq/data/encode/process";
my $raw_stat_file = "encode_raw_stat.csv";

my $op = "insert_bed";
my $dryrun;

# insert_bed
my $style = "center_intact";
my $noclean;    # do not clean ofg tables

# filter by filename regex
my $regex;

# filter by item counts
my $count;

# filter by meta info
my ( @dataType, @cell, @antibody );

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
    'regex=s'      => \$regex,
    'count=i'      => \$count,
    'dataType=s'   => \@dataType,
    'cell=s'       => \@cell,
    'antibody=s'   => \@antibody,
    'op=s'         => \$op,
    'dryrun'       => \$dryrun,
    'style=s'      => \$style,
    'parallel=i'   => \$parallel,
    'batch=i'      => \$batch_number,
    'noclean'      => \$noclean,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# bed info
#----------------------------------------------------------#
$stopwatch->start_message("Processing file ");

my @ymls;
{
    my $csv = Text::CSV_XS->new( { binary => 1 } );
    $csv->eol("\n");
    open my $fh, "<", $raw_stat_file;
    while ( my $row = $csv->getline($fh) ) {
        my $yml = {};
        $yml->{dataType}     = $row->[0];
        $yml->{cell_tag}     = $row->[2];
        $yml->{antibody_tag} = $row->[4];
        $yml->{itemCount}    = $row->[5];
        $yml->{filename}     = $row->[7];
        push @ymls, $yml;
    }
    close $fh;
}

if ($regex) {
    @ymls = grep { $_->{filename} =~ /$regex/i } @ymls;
}
if ($count) {
    @ymls = grep { $_->{itemCount} >= $count } @ymls;
}
if (@dataType) {
    my $set = Set::Scalar->new;
    $set->insert($_) for @dataType;
    @ymls = grep { $set->has( $_->{dataType} ) } @ymls;
}
if (@cell) {
    my $set = Set::Scalar->new;
    $set->insert($_) for @cell;
    @ymls = grep { $set->has( $_->{cell_tag} ) } @ymls;
}
if (@antibody) {
    my $set = Set::Scalar->new;
    $set->insert($_) for @antibody;
    @ymls = grep { $set->has( $_->{antibody_tag} ) } @ymls;
}

printf "Filtered %d\n", scalar @ymls;

if ( $op eq "insert_bed" ) {
    my $cmd
        = "perl $FindBin::Bin/../ofg/insert_bed.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db"
        . " --style $style"
        . " --batch $batch_number"
        . " --parallel $parallel"
        . ( $noclean ? " --noclean" : "" );

    for my $yml (@ymls) {
        if ( $yml->{dataType} eq "DnaseSeq" or $yml->{dataType} eq "FaireSeq" )
        {
            $cmd
                .= " --tag " . $yml->{cell_tag} . " --file " . $yml->{filename};
        }
        elsif ( $yml->{dataType} eq "TFBS" or $yml->{dataType} eq "Histone" ) {
            $cmd
                .= " --tag "
                . $yml->{cell_tag}
                . " --type "
                . $yml->{antibody_tag}
                . " --file "
                . $yml->{filename};
        }
    }
    exec_cmd( $cmd, $dryrun );
}
elsif ( $op eq "merge_to_runlist" ) {
    my $name = "merge";
    if (@dataType) {
        $name .= "_" . join "_", uniq(@dataType);
    }
    if (@cell) {
        $name .= "_" . join "_", uniq(@cell);
    }
    if (@antibody) {
        $name .= "_" . join "_", uniq(@antibody);
    }
    $name .= ".yml";

    my $cmd
        = "perl $FindBin::Bin/../ofg/bed_op.pl" . " --op $op" . " --name $name";

    for my $yml (@ymls) {
        $cmd .= " --file " . $yml->{filename};
    }
    exec_cmd( $cmd, $dryrun );
}

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub exec_cmd {
    my $cmd = shift;

    print "\n", "-" x 12, "CMD", "-" x 15, "\n";
    print $cmd , "\n";
    print "-" x 30, "\n";

    system $cmd unless $dryrun;
}

__END__

=head1 SYNOPSIS

    batch_encode.pl [options]
      Options:
        --help            brief help message
        --man             full documentation
        --file            csv file name
        --rep             times of bootstrap simulations

=cut
