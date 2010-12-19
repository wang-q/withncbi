#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Text::CSV_XS;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

# Database init values
my $rep_file = $Config->{rep}->{file};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?' => \$help,
    'man'    => \$man,
    'file=s' => \$rep_file,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Generate rep files...");

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my ( %early, %late );

my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
open my $csv_fh, "<", $rep_file or die "$rep_file: $!";
$csv->getline($csv_fh);    # bypass title line
while ( my $row = $csv->getline($csv_fh) ) {
    my $chr   = $row->[0];
    my $begin = $row->[1];
    my $end   = $row->[2];

    my $one_sixth = int( ( $end - $begin ) / 6 );

    if ( exists $early{$chr} ) {
        $early{$chr}->add( "$begin-" . ( $begin + $one_sixth ) );
        $early{$chr}->add( ( $end - $one_sixth ) . "-$end" );
    }
    else {
        $early{$chr} = AlignDB::IntSpan->new;
        $early{$chr}->add( "$begin-" . ( $begin + $one_sixth ) );
        $early{$chr}->add( ( $end - $one_sixth ) . "-$end" );
    }

    if ( exists $late{$chr} ) {
        $late{$chr}->add(
            ( $begin + 2 * $one_sixth ) . "-" . ( $begin + 4 * $one_sixth ) );
    }
    else {
        $late{$chr} = AlignDB::IntSpan->new;
        $late{$chr}->add(
            ( $begin + 2 * $one_sixth ) . "-" . ( $begin + 4 * $one_sixth ) );
    }
}
close $csv_fh;

unless ( -e "early" ) {
    mkdir "early", 0777
        or die "Cannot create directory: $!";
}
for my $key ( sort keys %early ) {
    DumpFile( "early/$key.yaml", $early{$key}->runlist );
    print "$key: ", $early{$key}->size, "\n";
}

unless ( -e "late" ) {
    mkdir "late", 0777
        or die "Cannot create directory: $!";
}
for my $key ( sort keys %late ) {
    DumpFile( "late/$key.yaml", $late{$key}->runlist );
    print "$key: ", $late{$key}->size, "\n";
}

$stopwatch->end_message;
exit;

__END__

Genome Res.  2007.   17:  1278-1285 

Human gene organization driven by the coordination of replication and
transcription

http://genome.cshlp.org/content/17/9/1278.full
