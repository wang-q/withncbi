#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use IO::Zlib;
use Set::Scalar;
use Tie::IxHash;
use AlignDB::IntSpan;
use App::RL::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

rmout2runlist.pl - Convert RepeatMasker .out file to chromosome runlists

=head1 SYNOPSIS

    perl rmout2runlist.pl [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     one gff3 file, gzipped file is supported
        --size          -s  STR     chr.sizes

=cut

GetOptions(
    'help|?'   => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s' => \my $infile,
    'size|s=s' => \my $size,
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Processing
#----------------------------------------------------------#
# only keep chromosomes listed in chr.sizes
my $chr_set_of = {};
if ( defined $size ) {
    if ( path($size)->is_file ) {
        printf "==> Load chr.sizes\n";
        my $length_of = App::RL::Common::read_sizes($size);
        for my $chr ( keys %{$length_of} ) {
            my $set = AlignDB::IntSpan->new;
            $set->add_pair( 1, $length_of->{$chr} );
            $chr_set_of->{$chr} = $set;
        }
    }
    else {
        die "*** [$size] doesn't exist!\n";
    }
}

printf "==> Load file\n";
my $in_fh = IO::Zlib->new( $infile, "rb" );

# runlists
my $all_repeat = {};    # all repeats combined

tie my %qr_of, 'Tie::IxHash';
%qr_of = (
    DNA            => qr{^(DNA)\??\/?([\w-]*)\??$}i,
    LINE           => qr{^(LINE)\??\/?([\w-]*)\??$}i,
    LTR            => qr{^(LTR)\??\/?([\w-]*)\??$}i,
    SINE           => qr{^(SINE)\??\/?([\w-]*)\??$}i,
    RC             => qr{^(RC)\??\/?([\w-]*)\??$}i,
    Other          => qr{^(Other)\??\/?([\w-]*)\??$}i,
    Low_complexity => qr{^(Low_complexity)\??\/?([\w-]*)\??$}i,
    Satellite      => qr{^(Satellite)\??\/?([\w-]*)\??$}i,
    Simple_repeat  => qr{^(Simple_repeat)\??\/?([\w-]*)\??$}i,
    Unknown        => qr{^(Unknown)$}i,
    Small_RNA      => qr{^\w*(RNA)$}i,
    UNCLASSIFIED   => qr{^(UNCLASSIFIED)$}i,                      # not exists
);

while (1) {
    my $line = <$in_fh>;
    last unless $line;
    next unless $line =~ /^\s*\d+/;

    chomp $line;
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    my @array = split /\s+/, $line;
    unless ( @array == 15 or @array == 16 ) {
        print "Fields error\n    $line\n";
        next;
    }
    my $type = $array[10];

    my $chr   = $array[4];
    my $start = $array[5];
    my $end   = $array[6];

    # only keep chromosomes listed in chr.sizes
    if ( defined $size ) {
        next unless exists $chr_set_of->{$chr};
    }

    my ( $family, $subfamily );
    for my $c ( keys %qr_of ) {
        if ( $type =~ $qr_of{$c} ) {
            $family    = $c;
            $subfamily = $2;
            last;
        }
    }

    if ( !defined $family ) {
        print "$type\n";
        $family = "UNCLASSIFIED";
    }

    # initialize sets
    if ( !exists $all_repeat->{$family} ) {
        $all_repeat->{$family} = {};
        $all_repeat->{$family} = {};
    }
    if ( !exists $all_repeat->{$family}{$chr} ) {
        $all_repeat->{$family}{$chr} = AlignDB::IntSpan->new;
    }

    # add sets
    $all_repeat->{$family}{$chr}->add_pair( $start, $end );
}
$in_fh->close;

#----------------------------------------------------------#
# Outputs
#----------------------------------------------------------#
print "==> Write Output files\n";
for my $f ( keys %{$all_repeat} ) {
    for my $chr ( keys %{ $all_repeat->{$f} } ) {
        $all_repeat->{$f}{$chr} = $all_repeat->{$f}{$chr}->runlist;
    }
    DumpFile( "$f.yml", $all_repeat->{$f} );
}

exit;
