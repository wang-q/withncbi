#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use File::Basename;
use String::Compare;
use List::MoreUtils qw(zip);
use MCE::Flow;
use Set::Scalar;

use App::Fasops::Common;
use AlignDB::IntSpan;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new;

=head1 NAME

masked_chr.pl - Soft-masking fa file.

=head1 SYNOPSIS

    perl masked_chr.pl [options]
      Options:
        --help      -?          brief help message
        --file      -f  STR     use a yaml file as annotation
        --dir       -d  STR     .fa dir

    mkdir example
    cd example
    wget -N ftp://ftp.ensembl.org/pub/release-82/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

    faops split-name Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz genome
    faops size genome/*.fa > genome/chr.sizes

    perl ~/Scripts/alignDB/slice/write_runlist_feature.pl -e yeast --feature repeat

    # Or use gff2runlist.pl or rmout2runlist.pl
    # wget -N ftp://ftp.ensembl.org/pub/release-82/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.82.gff3.gz

    perl ~/Scripts/withncbi/taxon/masked_chr.pl --dir genome --file yeast.repeat.yml

=cut

my $linelen = 60;

GetOptions(
    'help|?'   => sub { Getopt::Long::HelpMessage(0) },
    'file|f=s' => \my $yaml_file,
    'dir|d=s'  => \my $dir_fa,
    'parallel=i' => \( my $parallel = 1 ),
) or Getopt::Long::HelpMessage(1);

#----------------------------------------------------------#
# Init objects
#----------------------------------------------------------#
$stopwatch->start_message("Write masked chr...");

my $ftr_of = YAML::Syck::LoadFile($yaml_file);

#----------------------------#
# Soft mask
#----------------------------#
my @files = sort File::Find::Rule->file->name('*.fa')->in($dir_fa);

my @chrs = map { basename $_ , ".fa" } @files;    # strip dir and suffix
while (1) {
    my $lcss = lcss(@chrs);
    last unless $lcss;
    print "LCSS [$lcss]\n";
    my $rx = quotemeta $lcss;
    $chrs[$_] =~ s/$rx// for 0 .. $#chrs;
}
my $file_of = { zip( @chrs, @files ) };

for my $ftr_chr ( keys %{$ftr_of} ) { }

my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $ftr_chr = $chunk_ref->[0];

    my $ftr_chr_cmp = $ftr_chr;
    if ( $ftr_chr_cmp =~ /^chr/ ) {
        $ftr_chr_cmp =~ s/chr//;
    }

    # use the most similar file name
    my ($file_chr) = map { $_->[0] }
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, compare( $_, $ftr_chr_cmp ) ] } keys %{$file_of};

    printf "Write masked seq for ftr_chr [%s]\tfile_chr [%s]\n", $ftr_chr, $file_chr;

    # feature set
    my $ftr_set = AlignDB::IntSpan->new( $ftr_of->{$ftr_chr} );

    # seq
    my $seq_of = App::Fasops::Common::read_fasta( $file_of->{$file_chr} );
    my $seq = $seq_of->{ ( keys %{$seq_of} )[0] };

    my @sets = $ftr_set->sets;
    for my $set (@sets) {
        my $offset = $set->min - 1;
        my $length = $set->size;

        my $str = substr $seq, $offset, $length;
        $str = lc $str;
        substr $seq, $offset, $length, $str;
    }

    open my $out_fh, '>', $file_of->{$file_chr} . ".masked";
    print {$out_fh} ">$ftr_chr\n";
    print {$out_fh} substr( $seq, 0, $linelen, '' ) . "\n" while ($seq);
    close $out_fh;

    MCE->gather( $file_of->{$file_chr} );
};

MCE::Flow::init {
    chunk_size  => 1,
    max_workers => $parallel,
};
my @matched_files = mce_flow $worker, [ sort keys %{$ftr_of} ];
MCE::Flow::finish;

{    # combine all unmatched filess to chrUn.fasta
    my $fa_set = Set::Scalar->new;
    $fa_set->insert($_) for @files;
    $fa_set->delete($_) for @matched_files;

    print "\n";
    printf "There are %d unmatched file(s)\n", $fa_set->size;
    if ( $fa_set->size ) {
        my $str = join " ", map { basename $_} $fa_set->elements;
        print "We'll combine these files into chrUn.fasta\n";
        print $str, "\n";

        system "cat $str > chrUn.fasta";
        system "rm $str";
    }
}

$stopwatch->end_message;

# comes from
# http://stackoverflow.com/questions/499967/how-do-i-determine-the-longest-similar-portion-of-several-strings
sub lcss {
    return '' unless @_;
    return $_[0] if @_ == 1;
    my $i          = 0;
    my $first      = shift;
    my $min_length = length($first);
    for (@_) {
        $min_length = length($_) if length($_) < $min_length;
    }
INDEX: for my $ch ( split //, $first ) {
        last INDEX unless $i < $min_length;
        for my $string (@_) {
            last INDEX if substr( $string, $i, 1 ) ne $ch;
        }
    }
    continue { $i++ }
    return substr $first, 0, $i;
}

__END__
