#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use IO::Zlib;
use AlignDB::IntSpan;

use lib "$FindBin::RealBin/../lib";
use MyUtil qw(read_sizes);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

gff2runlist.pl - Convert gff3 file to chromosome runlists

=head1 SYNOPSIS

    perl gff2runlist.pl [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     one gff3 file, gzipped file is supported
        --range         -r  INT     range of up/down-stream. Default is [5000]
        --size          -s  STR     chr.sizes

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \my $infile,
    'range|r=i' => \( my $range = 5000 ),
    'size|s=s' => \my $size,
) or HelpMessage(1);

#----------------------------------------------------------#
# Processing
#----------------------------------------------------------#
# avoid out of chromosome for up/down-stream
my $chr_set_of = {};
if ( defined $size and path($size)->is_file ) {
    printf "Load chr.sizes\n";
    my $length_of = read_sizes($size);
    for my $chr ( keys %{$length_of} ) {
        my $set = AlignDB::IntSpan->new;
        $set->add_pair( 1, $length_of->{$chr} );
        $chr_set_of->{$chr} = $set;
    }
}

my $in_fh = IO::Zlib->new( $infile, "rb" );

# runlists
my $gene       = {};
my $upstream   = {};
my $downstream = {};
my $transcript = {};
my $intron     = {};
my $fiveutr    = {};
my $threeutr   = {};
my $cds        = {};

# current transcript
my $cur_transcript;
my @exon;
my @five;
my @three;
my @cds;

while ( my $line = <$in_fh> ) {
    next if $line =~ /^#/;
    chomp $line;
    my @array = split( "\t", $line );
    my $type = $array[2];

    my $chr    = $array[0];
    my $start  = $array[3];
    my $end    = $array[4];
    my $strand = $array[6];    # strand may be "."
    if ( $strand ne "-" ) {
        $strand = "+";
    }

    my @attrs = split ";", $array[8];
    my %attr_of;
    for my $item (@attrs) {
        my @pair = split "=", $item;
        if ( @pair == 2 ) {
            $attr_of{ $pair[0] } = $pair[1];
        }
    }

    if ( $type eq 'gene' ) {
        my $gene_id;
        if ( exists $attr_of{gene_id} ) {
            $gene_id = $attr_of{gene_id};
        }
        else {
            $gene_id = $attr_of{ID};
            print " " x 4 . "Can't get gene_id\n";
            print " " x 4 . "$line\n";
        }

        # gene
        printf "gene: %s\n", $gene_id;
        my $set = AlignDB::IntSpan->new;
        $set->add_pair( $start, $end );
        $gene->{$gene_id} = { $chr => $set->runlist };

        # up/down-stream
        my $up_set = AlignDB::IntSpan->new;
        $up_set->add_pair( $start - $range, $start - 1 );

        my $down_set = AlignDB::IntSpan->new;
        $down_set->add_pair( $end + 1, $end + $range );

        if ( $strand eq "-" ) {
            ( $up_set, $down_set ) = ( $down_set, $up_set );
        }
        if ( exists $chr_set_of->{$chr} ) {
            $up_set   = $up_set->intersect( $chr_set_of->{$chr} );
            $down_set = $down_set->intersect( $chr_set_of->{$chr} );
        }
        $upstream->{$gene_id}   = { $chr => $up_set->runlist };
        $downstream->{$gene_id} = { $chr => $down_set->runlist };
    }

    if ( ( defined $cur_transcript ) and ( $type eq "transcript" ) ) {
        printf "transcript: %s\n", $cur_transcript;

        # process collected features
        my $exon_set = AlignDB::IntSpan->new;
        for (@exon) {
            $exon_set->add_pair( $_->[1], $_->[2] );
        }

        my $five_set = AlignDB::IntSpan->new;
        for (@five) {
            $five_set->add_pair( $_->[1], $_->[2] );
        }

        my $three_set = AlignDB::IntSpan->new;
        for (@three) {
            $three_set->add_pair( $_->[1], $_->[2] );
        }

        my $cds_set = AlignDB::IntSpan->new;
        for (@cds) {
            $cds_set->add_pair( $_->[1], $_->[2] );
        }

        $transcript->{$cur_transcript} = { $chr => $exon_set->runlist };
        $intron->{$cur_transcript}     = { $chr => $exon_set->holes->runlist };
        $fiveutr->{$cur_transcript}    = { $chr => $five_set->runlist };
        $threeutr->{$cur_transcript}   = { $chr => $three_set->runlist };
        $cds->{$cur_transcript}        = { $chr => $cds_set->runlist };

        # initialize the next transcript
        my $transcript_id;
        if ( exists $attr_of{transcript_id} ) {
            $transcript_id = $attr_of{transcript_id};
        }
        else {
            $transcript_id = $attr_of{ID};
            print " " x 4 . "Can't get transcript_id\n";
            print " " x 4 . "$line\n";
        }
        $cur_transcript = $transcript_id;

        # empty caches
        @exon  = ();
        @five  = ();
        @three = ();
        @cds   = ();
    }
    elsif ( $type eq "transcript" ) {    # First transcript
        my $transcript_id;
        if ( exists $attr_of{transcript_id} ) {
            $transcript_id = $attr_of{transcript_id};
        }
        else {
            $transcript_id = $attr_of{ID};
            print " " x 4 . "Can't get transcript_id\n";
            print " " x 4 . "$line\n";
        }
        $cur_transcript = $transcript_id;
    }
    elsif ( $type eq 'exon' ) {
        push @exon, [ $chr, $start, $end ];
    }
    elsif ( $type eq 'five_prime_UTR' ) {
        push @five, [ $chr, $start, $end ];
    }
    elsif ( $type eq 'three_prime_UTR' ) {
        push @three, [ $chr, $start, $end ];
    }
    elsif ( $type eq 'CDS' ) {
        push @cds, [ $chr, $start, $end ];
    }
}

# last transcript
if ( scalar @exon > 0 ) {
    printf "transcript: %s\n", $cur_transcript;

    # process collected features
    my $exon_set = AlignDB::IntSpan->new;
    for (@exon) {
        $exon_set->add_pair( $_->[1], $_->[2] );
    }

    my $five_set = AlignDB::IntSpan->new;
    for (@five) {
        $five_set->add_pair( $_->[1], $_->[2] );
    }

    my $three_set = AlignDB::IntSpan->new;
    for (@three) {
        $three_set->add_pair( $_->[1], $_->[2] );
    }

    my $cds_set = AlignDB::IntSpan->new;
    for (@cds) {
        $cds_set->add_pair( $_->[1], $_->[2] );
    }

    $transcript->{$cur_transcript} = { $exon[0]->[0] => $exon_set->runlist };
    $intron->{$cur_transcript}     = { $exon[0]->[0] => $exon_set->holes->runlist };
    $fiveutr->{$cur_transcript}    = { $exon[0]->[0] => $five_set->runlist };
    $threeutr->{$cur_transcript}   = { $exon[0]->[0] => $three_set->runlist };
    $cds->{$cur_transcript}        = { $exon[0]->[0] => $cds_set->runlist };
}

$in_fh->close;

DumpFile( "gene.yml",       $gene );
DumpFile( "upstream.yml",   $upstream );
DumpFile( "downstream.yml", $downstream );
DumpFile( "transcript.yml", $transcript );
DumpFile( "intron.yml",     $intron );
DumpFile( "fiveutr.yml",    $fiveutr );
DumpFile( "threeutr.yml",   $threeutr );
DumpFile( "cds.yml",        $cds );

exit;
