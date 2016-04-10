#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Path::Tiny;
use IO::Zlib;
use Set::Scalar;
use Tie::IxHash;
use AlignDB::IntSpan;

use lib "$FindBin::RealBin/../lib";
use MyUtil qw(read_sizes);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

rmout2runlist.pl - Convert RepeatMasker .out file to chromosome runlists

=head1 SYNOPSIS

    perl gff2runlist.pl [options]
      Options:
        --help          -?          brief help message
        --file          -f  STR     one gff3 file, gzipped file is supported
        --size          -s  STR     chr.sizes

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \my $infile,
    'size|s=s' => \my $size,
) or HelpMessage(1);

#----------------------------------------------------------#
# Processing
#----------------------------------------------------------#
# only keep chromosomes listed in chr.sizes
my $chr_set_of = {};
if ( defined $size ) {
    if ( path($size)->is_file ) {
        printf "==> Load chr.sizes\n";
        my $length_of = read_sizes($size);
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
my $all_gene = {};    # all genes combined

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

my %count_of = map { $_ => undef } keys %qr_of;

# chromosome names
my $all_name_set = Set::Scalar->new;

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
            $count_of{$family}++;
            last;
        }
    }

    if ( !defined $family ) {
        print "$type\n";
        $count_of{UNCLASSIFIED}++;
    }
}
$in_fh->close;

print Dump \%count_of;

# #----------------------------------------------------------#
# # Outputs
# #----------------------------------------------------------#
# print "==> Write Output files\n";
# for my $f (qw{gene upstream downstream exon five_prime_UTR three_prime_UTR CDS intron}) {
#     for my $chr ( $all_name_set->members ) {
#         $all_gene->{$f}{$chr} = $all_gene->{$f}{$chr}->runlist;
#     }
#     for my $id ( keys %{ $sep_gene->{$f} } ) {
#         for my $chr ( keys %{ $sep_gene->{$f}{$id} } ) {
#             $sep_gene->{$f}{$id}{$chr} = $sep_gene->{$f}{$id}{$chr}->runlist;
#         }
#     }
#     DumpFile( "sep-$f.yml", $sep_gene->{$f} );
# }
# DumpFile( "all-gene.yml", $all_gene );

exit;
