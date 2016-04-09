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
        --size          -s  STR     chr.sizes
        --range         -r  INT     range of up/down-stream. Default is [1000]
        --clean         -c          up/down-stream don't contain any regions of other genes

=cut

GetOptions(
    'help|?'   => sub { HelpMessage(0) },
    'file|f=s' => \my $infile,
    'size|s=s' => \my $size,
    'range|r=i' => \( my $range = 1000 ),
    'clean|c' => \my $clean,
) or HelpMessage(1);

#----------------------------------------------------------#
# Processing
#----------------------------------------------------------#
# only keep chromosomes listed in chr.sizes
# avoid out of chromosome for up/down-stream
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
my $sep_gene = {};    # seperated genes

# chromosome names
my $all_name_set = Set::Scalar->new;

# current transcript
my $cur_transcript;
my $cache = {
    exon            => [],
    five_prime_UTR  => [],
    three_prime_UTR => [],
    CDS             => [],
};

while (1) {
    my $line = <$in_fh>;
    last unless $line;
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

    # only keep chromosomes listed in chr.sizes
    if ( defined $size ) {
        next unless exists $chr_set_of->{$chr};
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
        printf "gene: %s\n", $gene_id;

        $all_name_set->insert($chr);

        # gene and up/down-stream
        my $set_of = {
            gene       => AlignDB::IntSpan->new,
            upstream   => AlignDB::IntSpan->new,
            downstream => AlignDB::IntSpan->new,
        };
        $set_of->{gene}->add_pair( $start, $end );
        $set_of->{upstream}->add_pair( $start - $range, $start - 1 );
        $set_of->{downstream}->add_pair( $end + 1, $end + $range );

        # swap up/down-stream for negtive strand
        if ( $strand eq "-" ) {
            ( $set_of->{upstream}, $set_of->{downstream} )
                = ( $set_of->{downstream}, $set_of->{upstream} );
        }

        # avoid out of chromosome for up/down-stream
        if ( exists $chr_set_of->{$chr} ) {
            $set_of->{upstream}   = $set_of->{upstream}->intersect( $chr_set_of->{$chr} );
            $set_of->{downstream} = $set_of->{downstream}->intersect( $chr_set_of->{$chr} );
        }

        for my $f (qw{gene upstream downstream}) {

            # initialize sets
            if ( !exists $all_gene->{$f} ) {
                $all_gene->{$f} = {};
                $sep_gene->{$f} = {};
            }
            if ( !exists $all_gene->{$f}{$chr} ) {
                $all_gene->{$f}{$chr} = AlignDB::IntSpan->new;
            }

            # add sets
            $all_gene->{$f}{$chr}->add( $set_of->{$f} );
            $sep_gene->{$f}{$gene_id} = { $chr => $set_of->{$f} };
        }
    }

    if ( ( defined $cur_transcript ) and ( $type eq "transcript" ) ) {
        printf "transcript: %s\n", $cur_transcript;

        process_transcript( $all_gene, $sep_gene, $cur_transcript, $cache );

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
        $cache = {
            exon            => [],
            five_prime_UTR  => [],
            three_prime_UTR => [],
            CDS             => [],
        };
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
        push @{ $cache->{$type} }, [ $chr, $start, $end ];
    }
    elsif ( $type eq 'five_prime_UTR' ) {
        push @{ $cache->{$type} }, [ $chr, $start, $end ];
    }
    elsif ( $type eq 'three_prime_UTR' ) {
        push @{ $cache->{$type} }, [ $chr, $start, $end ];
    }
    elsif ( $type eq 'CDS' ) {
        push @{ $cache->{$type} }, [ $chr, $start, $end ];
    }
}
$in_fh->close;

# last transcript
if ( scalar @{ $cache->{exon} } > 0 ) {
    printf "transcript: %s\n", $cur_transcript;
    process_transcript( $all_gene, $sep_gene, $cur_transcript, $cache );
}

# intergenic regions
{
    print "==> Intergenic\n";
    $all_gene->{intergenic} = {};
    for my $chr ( $all_name_set->members ) {
        my $set = $chr_set_of->{$chr}->copy;
        for my $f (qw{gene upstream downstream}) {
            if ( exists $all_gene->{$f}{$chr} ) {
                $set->remove( $all_gene->{$f}{$chr} );
            }
        }
        $all_gene->{intergenic}{$chr} = $set->runlist;
    }
}

# clean up/down-stream
if ($clean) {
    for my $f (qw{upstream downstream}) {
        print "==> Clean $f\n";
        for my $chr ( $all_name_set->members ) {
            $all_gene->{$f}{$chr}->remove( $all_gene->{gene}{$chr} );
        }
        for my $id ( keys %{ $sep_gene->{gene} } ) {
            for my $chr ( keys %{ $sep_gene->{gene}{$id} } ) {

                #print "$f $id $chr\n";
                $sep_gene->{$f}{$id}{$chr}->remove( $all_gene->{gene}{$chr} );
            }
        }
    }
}

#----------------------------------------------------------#
# Outputs
#----------------------------------------------------------#
print "==> Write Output files\n";
for my $f (qw{gene upstream downstream exon five_prime_UTR three_prime_UTR CDS intron}) {
    for my $chr ( $all_name_set->members ) {
        $all_gene->{$f}{$chr} = $all_gene->{$f}{$chr}->runlist;
    }
    for my $id ( keys %{ $sep_gene->{$f} } ) {
        for my $chr ( keys %{ $sep_gene->{$f}{$id} } ) {
            $sep_gene->{$f}{$id}{$chr} = $sep_gene->{$f}{$id}{$chr}->runlist;
        }
    }
    DumpFile( "sep-$f.yml", $sep_gene->{$f} );
}
DumpFile( "all-gene.yml", $all_gene );

exit;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#
sub process_transcript {
    my $all_gene       = shift;
    my $sep_gene       = shift;
    my $cur_transcript = shift;
    my $cache          = shift;

    # process cached features
    my $set_of = {
        exon            => AlignDB::IntSpan->new,
        five_prime_UTR  => AlignDB::IntSpan->new,
        three_prime_UTR => AlignDB::IntSpan->new,
        CDS             => AlignDB::IntSpan->new,
        intron          => AlignDB::IntSpan->new,
    };

    for my $f (qw{exon five_prime_UTR three_prime_UTR CDS}) {
        for my $item ( @{ $cache->{$f} } ) {
            $set_of->{$f}->add_pair( $item->[1], $item->[2] );
        }
    }
    $set_of->{intron} = $set_of->{exon}->holes;

    my $chr = $cache->{exon}[0][0];
    for my $f (qw{exon five_prime_UTR three_prime_UTR CDS intron}) {

        # initialize sets
        if ( !exists $all_gene->{$f} ) {
            $all_gene->{$f} = {};
            $sep_gene->{$f} = {};
        }
        if ( !exists $all_gene->{$f}{$chr} ) {
            $all_gene->{$f}{$chr} = AlignDB::IntSpan->new;
        }

        # add sets
        $all_gene->{$f}{$chr}->add( $set_of->{$f} );
        $sep_gene->{$f}{$cur_transcript} = { $chr => $set_of->{$f} };
    }

    return;
}
