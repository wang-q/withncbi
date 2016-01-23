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
use File::Spec;

use Text::CSV_XS;

use AlignDB::Run;
use AlignDB::Stopwatch;

use FindBin;

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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

# executable file location
my $info_file = "/home/wangq/data/encode/stat/encode_raw_stat.csv";

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
my ( @dataType, @cell, @antibody, @view );

my $filename;

my $remove_chr;

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
    'i|info=s'     => \$info_file,
    'regex=s'      => \$regex,
    'count=i'      => \$count,
    'dataType=s'   => \@dataType,
    'cell=s'       => \@cell,
    'antibody=s'   => \@antibody,
    'view=s'       => \@view,
    'op=s'         => \$op,
    'dryrun'       => \$dryrun,
    'n|name=s'     => \$filename,
    'r|remove'     => \$remove_chr,
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
    open my $fh, "<", $info_file;
    while ( my $row = $csv->getline($fh) ) {
        my $yml = {};
        $yml->{dataType}     = $row->[0];
        $yml->{cell}         = $row->[1];
        $yml->{cell_tag}     = $row->[2];
        $yml->{antibody}     = $row->[3];
        $yml->{antibody_tag} = $row->[4];
        $yml->{view}         = $row->[5];
        $yml->{itemCount}    = $row->[6];
        $yml->{average_size} = $row->[7];
        $yml->{filename}     = $row->[8];
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
if (@view) {
    my $set = Set::Scalar->new;
    $set->insert($_) for @view;
    @ymls = grep { $set->has( $_->{view} ) } @ymls;
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
        if (   $yml->{dataType} eq "DnaseSeq"
            or $yml->{dataType} eq "FaireSeq"
            or $yml->{dataType} eq "RepliChip" )
        {
            $cmd
                .= "\\\n    --tag " . $yml->{cell_tag} . " --file " . $yml->{filename};
        }
        elsif ( $yml->{dataType} eq "TFBS" or $yml->{dataType} eq "Histone" ) {
            $cmd
                .= "\\\n    --tag "
                . $yml->{cell_tag}
                . " --type "
                . $yml->{antibody_tag}
                . " --file "
                . $yml->{filename};
        }
        elsif ( $yml->{dataType} eq "RepliSeq" ) {
            $cmd
                .= "\\\n    --tag "
                . $yml->{cell_tag}
                . " --type "
                . $yml->{view}
                . " --file "
                . $yml->{filename};
        }
    }
    exec_cmd( $cmd, $dryrun );
}
elsif ( $op eq "count_bed" ) {
    my $cmd
        = "perl $FindBin::Bin/../ofg/count_bed.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db"
        . " --batch $batch_number"
        . " --parallel $parallel"
        . ( $noclean ? " --noclean" : "" );

    for my $yml (@ymls) {
        if (   $yml->{dataType} eq "DnaseSeq"
            or $yml->{dataType} eq "FaireSeq"
            or $yml->{dataType} eq "RepliChip" )
        {
            $cmd .= " --file " . $yml->{filename};
        }
        elsif ( $yml->{dataType} eq "TFBS" or $yml->{dataType} eq "Histone" ) {
            $cmd .= " --file " . $yml->{filename};
        }
    }
    exec_cmd( $cmd, $dryrun );
}
elsif ( $op eq "merge_to_runlist" ) {
    if ( !$filename ) {
        $filename = "merge";
        if (@dataType) {
            $filename .= "_" . join "_", uniq(@dataType);
        }
        if (@cell) {
            $filename .= "_" . join "_", uniq(@cell);
        }
        if (@antibody) {
            $filename .= "_" . join "_", uniq(@antibody);
        }
        $filename .= ".yml";
    }

    my $cmd
        = "perl $FindBin::Bin/../ofg/bed_op.pl"
        . " --op $op"
        . " --name $filename";
    $cmd .= " --remove " if $remove_chr;

    for my $yml (@ymls) {
        $cmd .= " --file " . $yml->{filename};
    }
    exec_cmd( $cmd, $dryrun );
}
elsif ( $op eq "bed_diff" or $op eq "bed_intersect" ) {
    my $this_info_file = $op;
    if (@dataType) {
        $this_info_file .= "_" . join "_", uniq(@dataType);
    }
    if (@cell) {
        $this_info_file .= "_" . join "_", uniq(@cell);
    }
    if (@antibody) {
        $this_info_file .= "_" . join "_", uniq(@antibody);
    }
    $this_info_file .= ".csv";
    open my $fh, '>', $this_info_file;
    for my $yml (@ymls) {
        my $newname = File::Spec->catfile( File::Spec->curdir,
            basename( $yml->{filename} ) . ".$op" );
        $newname = File::Spec->rel2abs($newname);
        print {$fh} $yml->{dataType} . ","
            . $yml->{cell} . ","
            . $yml->{cell_tag} . ","
            . $yml->{antibody} . ","
            . $yml->{antibody_tag} . ","
            . $yml->{itemCount} . ","
            . $yml->{average_size} . ","
            . $newname . "," . "\n";
    }
    close $fh;

    my $worker = sub {
        my $job = shift;
        my $cmd
            = "perl $FindBin::Bin/../ofg/bed_op.pl"
            . " --op $op"
            . " --name $filename";

        $cmd .= " --file " . $job->{filename};
        $cmd .= " --remove " if $remove_chr;
        exec_cmd( $cmd, $dryrun );
    };

    my $run = AlignDB::Run->new(
        parallel => $parallel,
        jobs     => \@ymls,
        code     => $worker,
    );
    $run->run;

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
