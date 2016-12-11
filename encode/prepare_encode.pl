#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use Digest::MD5;

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

# executable file location
my $kent_bin      = "/home/wangq/bin/x86_64";
my $file_metainfo = "/home/wangq/data/encode/files.txt";
my $dir_result    = "/home/wangq/data/encode/process";
my $file_chr_size = "/home/wangq/data/encode/human.chr.sizes";
my $file_encode;

my $filetype = "bigBed";    # bigWig, bed

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'f|file=s'   => \$file_encode,
    'm|meta=s'   => \$file_metainfo,
    'r|result=s' => \$dir_result,
    's|size=s'   => \$file_chr_size,
    't|type=s'   => \$filetype,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# bed info
#----------------------------------------------------------#
$stopwatch->start_message("Processing file [$file_encode]");

if ( !-e $file_encode ) {
    die "[$file_encode] does not exist.\n";
}

# make dirs
unless ( -e $dir_result ) {
    mkdir $dir_result, 0777;
}

my $basename = basename($file_encode);

if ( $filetype eq "bigBed" ) {

    # bigBedInfo - Show information about a bigBed file.
    # usage:
    #   bigBedInfo file.bb
    # options:
    #   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs
    #   -chroms - list all chromosomes and their sizes
    #   -zooms - list all zoom levels and theier sizes
    #   -as - get autoSql spec
    #
    # version: 4
    # isCompressed: yes
    # isSwapped: 0
    # itemCount: 167,357
    # primaryDataSize: 1,574,077
    # primaryIndexSize: 21,196
    # zoomLevels: 8
    # chromCount: 24
    # basesCovered: 25,094,170
    # meanDepth (of bases covered): 1.000374
    # minDepth: 1.000000
    # maxDepth: 2.000000
    # std of depth: 0.019330

    my $cmd = "$kent_bin/bigBedInfo" . " $file_encode" . " > $dir_result/$basename.info.txt";
    system $cmd if !-e "$dir_result/$basename.info.txt";

    # bigBedToBed - Convert from bigBed to ascii bed format.
    # usage:
    #   bigBedToBed input.bb output.bed
    # options:
    #   -chrom=chr1 - if set restrict output to given chromosome
    #   -start=N - if set, restrict output to only that over start
    #   -end=N - if set, restict output to only that under end
    #   -maxItems=N - if set, restrict output to first N items
    #   -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs
    $cmd = "$kent_bin/bigBedToBed" . " $file_encode " . " $dir_result/$basename.bed";
    system $cmd if !-e "$dir_result/$basename.bed";
}
elsif ( $filetype eq "bed" ) {
    print "    basename $basename\n";

    my $cmd;
    if ( $file_encode =~ /\.gz$/ ) {
        $cmd = "gunzip -c " . " $file_encode" . " > $dir_result/$basename.bed";
        system $cmd if !-e "$dir_result/$basename.bed";
    }
    else {
        $cmd = "cp " . " $file_encode" . " $dir_result/$basename.bed";
        system $cmd if !-e "$dir_result/$basename.bed";
    }

# bedToBigBed v. 4 - Convert bed file to bigBed.
# usage:
#    bedToBigBed in.bed chrom.sizes out.bb
# Where in.bed is in one of the ascii bed formats, but not including track lines
# and chrom.sizes is two column: <chromosome name> <size in bases>
# and out.bb is the output indexed big bed file.
# The in.bed file must be sorted by chromosome,start,
#   to sort a bed file, use the unix sort command:
#      sort -k1,1 -k2,2n unsorted.bed > sorted.bed
#
# options:
#    -blockSize=N - Number of items to bundle in r-tree.  Default 256
#    -itemsPerSlot=N - Number of data points bundled at lowest level. Default 512
#    -bedFields=N - Number of fields that fit standard bed definition.  If undefined
#                   assumes all fields in bed are defined.
#    -as=fields.as - If have non-standard fields, it's great to put a definition
#                    of each field in a row in AutoSql format here.
#    -unc - If set, do not use compression.   -tabs - If set, expect fields to be tab separated, normally
#            expects white space separator.
    $cmd
        = "$kent_bin/bedToBigBed"
        . " -bedFields=3"
        . " $dir_result/$basename.bed"
        . " $file_chr_size"
        . " $dir_result/$basename.bb";
    system $cmd if !-e "$dir_result/$basename.bb";

    $cmd
        = "$kent_bin/bigBedInfo"
        . " $dir_result/$basename.bb"
        . " > $dir_result/$basename.info.txt";
    system $cmd if !-e "$dir_result/$basename.info.txt";
}
elsif ( $filetype eq "bigWig" ) {

    # bigWigInfo - Print out information about bigWig file.
    # usage:
    #    bigWigInfo file.bw
    # options:
    #    -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs
    #    -chroms - list all chromosomes and their sizes
    #    -zooms - list all zoom levels and their sizes
    #    -minMax - list the min and max on a single line
    #
    # version: 4
    # isCompressed: yes
    # isSwapped: 0
    # primaryDataSize: 19,758,044
    # primaryIndexSize: 78,528
    # zoomLevels: 9
    # chromCount: 25
    # basesCovered: 129,721,134
    # mean: 0.047313
    # min: -1.765402
    # max: 1.942196
    # std: 0.853118
    my $cmd = "$kent_bin/bigWigInfo" . " $file_encode" . " > $dir_result/$basename.info.txt";
    system $cmd if !-e "$dir_result/$basename.info.txt";

    # bigWigToWig - Convert bigWig to wig.  This will keep more of the same structure of the
    # original wig than bigWigToBedGraph does, but still will break up large stepped sections
    # into smaller ones.
    # usage:
    #    bigWigToWig in.bigWig out.wig
    # options:
    #    -chrom=chr1 - if set restrict output to given chromosome
    #    -start=N - if set, restrict output to only that over start
    #    -end=N - if set, restict output to only that under end
    #    -udcDir=/dir/to/cache - place to put cache for remote bigBed/bigWigs
    $cmd = "$kent_bin/bigWigToWig" . " $file_encode " . " $dir_result/$basename.wig";
    system $cmd if !-e "$dir_result/$basename.wig";

    $cmd
        = "wig2bed --do-not-sort" . " < $dir_result/$basename.wig" . " > $dir_result/$basename.bed";
    system $cmd if !-e "$dir_result/$basename.bed";
}
else {
    die "Unknown filetype: $filetype\n";
}

if ( !-e "$dir_result/$basename.yml" ) {

    # read meta info
    open my $meta_fh, '<', $file_metainfo;
    my @lines;
    while ( my $line = <$meta_fh> ) {
        chomp $line;
        if ( index( $line, $basename ) != -1 ) {
            push @lines, $line;
        }
    }
    close $meta_fh;

    if ( scalar @lines == 0 ) {
        print "Can't find metainfo for $file_encode\n";
        exit;
    }
    elsif ( scalar @lines > 1 ) {
        print "Find multiple metainfo for $file_encode\n";
        print scalar @lines;
        exit;
    }
    elsif ( scalar @lines == 1 ) {
        print "Find metainfo for $file_encode\n";
    }

    my @fields = split /;\s+/, ( split /\t/, $lines[0] )[1];
    my %meta;
    for (@fields) {
        my ( $key, $value ) = split /=/;
        $meta{$key} = $value;
    }

    # add bedinfo
    open my $info_fh, '<', "$dir_result/$basename.info.txt";
    while ( my $line = <$info_fh> ) {
        chomp $line;
        if ( $line =~ /^(itemCount|basesCovered):\s+([\d,]+)$/ ) {
            my $key   = $1;
            my $value = $2;
            $value =~ s/,//g;
            $meta{$key} = $value;
        }
    }
    close $info_fh;

    if ( !exists $meta{itemCount} ) {
        ( $meta{itemCount} ) = split /\s+/, `wc -l $dir_result/$basename.bed`;
    }
    $meta{average_size} = $meta{basesCovered} / $meta{itemCount};

    # bed file
    $meta{filename} = "$dir_result/$basename.bed";

    DumpFile( "$dir_result/$basename.yml", \%meta );

    # MD5 checksum
    open my $bb_fh, '<', $file_encode;
    binmode($bb_fh);

    my $md5 = Digest::MD5->new->addfile($bb_fh)->hexdigest;
    close $bb_fh;

    if ( $md5 ne $meta{md5sum} ) {
        printf "Orginal:%s\t\nCalculated:%s\t\n", $meta{md5sum}, $md5;
        die "MD5 sum doesn't match\n";
    }
    else {
        print "MD5 sum OK\n";
    }
}

__END__

=head1 SYNOPSIS

    prepare_encode.pl [options]
      Options:
        --help            brief help message
        --man             full documentation
        --file            csv file name
        --rep             times of bootstrap simulations

=cut
