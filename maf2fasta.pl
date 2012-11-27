#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Spec;
use File::Find::Rule;
use File::Basename;
use IO::Zlib;

use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);
use AlignDB::Run;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_dir = '.';    # Specify location here
my $out_dir;         # Specify output dir here
my $length_threshold = 1000;    # Set the threshold of alignment length
my $has_outgroup;    # Treate last sequence as outgroup and move it to top
my $subset;          # get sequences of listed names, seperated by comma
my $block;           # write galaxy style blocked fasta

# run in parallel mode
my $parallel = 1;

my $gzip;            # open .maf.gz

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'i|in_dir=s'   => \$in_dir,
    'o|out_dir=s'  => \$out_dir,
    'length=i'     => \$length_threshold,
    'has_outgroup' => \$has_outgroup,
    'subset=s'     => \$subset,
    'block'        => \$block,
    'parallel=i'   => \$parallel,
    'gzip'         => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# make output dir
#----------------------------------------------------------#
if ($subset) {
    if ($has_outgroup) {
        print "--subset assigned, --has_outgroup will be ingored\n";
    }
    $has_outgroup = 0;
}

unless ($out_dir) {
    $out_dir = File::Spec->rel2abs($in_dir) . "_fasta";
    $out_dir = $out_dir . "_block" if $block;
    $out_dir = $out_dir . "_$subset" if $subset;
    $out_dir = $out_dir . "_ref_first" if $has_outgroup;
}
if ( !-e $out_dir ) {
    mkdir $out_dir, 0777;
}
elsif ($block) {
    print "We are going to output blocked fasta.\n";
    die "$out_dir exists, you should remove it first to avoid errors.\n";
}

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my @files;
if ( !$gzip ) {
    @files = sort File::Find::Rule->file->name('*.maf')->in($in_dir);
    printf "\n----Total .maf Files: %4s----\n\n", scalar @files;
}
if ( scalar @files == 0 or $gzip ) {
    @files = sort File::Find::Rule->file->name('*.maf.gz')->in($in_dir);
    printf "\n----Total .maf.gz Files: %4s----\n\n", scalar @files;
    $gzip++;
}

#----------------------------------------------------------#
# run
#----------------------------------------------------------#
my $worker = sub {
    my $infile = shift;

    my $in_fh;
    if ( !$gzip ) {
        open $in_fh, '<', $infile;
    }
    else {
        $in_fh = IO::Zlib->new( $infile, "rb" );
    }

    my $content = '';
ALN: while ( my $line = <$in_fh> ) {
        if ( $line =~ /^\s+$/ and $content =~ /\S/ ) {    # meet blank line
            my @slines = grep {/\S/} split /\n/, $content;
            $content = '';

            # parse maf
            my @names;
            my $info_of = {};
            for my $sline (@slines) {
                my ( $s, $src, $start, $size, $strand, $srcsize, $text )
                    = split /\s+/, $sline;

                my ( $species, $chr_name ) = split /\./, $src;
                $chr_name = $species if !defined $chr_name;

                # adjust coordinates to be one-based inclusive
                $start = $start + 1;

                push @names, $species;
                $info_of->{$species} = {
                    seq        => $text,
                    name       => $species,
                    chr_name   => $chr_name,
                    chr_start  => $start,
                    chr_end    => $start + $size - 1,
                    chr_strand => $strand,
                };
            }

            if ($subset) {
                my $name_str = join " ", @names;
                my @subsets = split ",", $subset;
                next ALN if scalar @names < scalar @subsets;
                for (@subsets) {
                    next ALN if $name_str !~ /$_/;
                }
                @names = @subsets;
            }
            my $target = $names[0];

            # if $subset, this will never been executed
            if ($has_outgroup) {
                my $outgoup = pop @names;
                unshift @names, $outgoup;
            }

            # output
            if ($block) {
                my $outfile = basename($infile);
                $outfile = $out_dir . "/$outfile" . ".fas";

                open my $fh, ">>", $outfile;
                for my $species (@names) {
                    my $chr    = $info_of->{$species}{chr_name};
                    my $start  = $info_of->{$species}{chr_start};
                    my $end    = $info_of->{$species}{chr_end};
                    my $strand = $info_of->{$species}{chr_strand};
                    print {$fh} ">$species.$chr($strand)";
                    print {$fh} ":$start-$end";
                    print {$fh} "|species=$species";
                    print {$fh} "\n";
                    print {$fh} $info_of->{$species}{seq}, "\n";
                }
                print {$fh} "\n";
                close $fh;
            }
            else {
                my $chr   = $info_of->{$target}{chr_name};
                my $start = $info_of->{$target}{chr_start};
                my $end   = $info_of->{$target}{chr_end};

                open my $fh, '>', $out_dir . "/$chr" . "-$start" . "-$end.fas";
                for my $species (@names) {
                    print {$fh} ">$species\n";
                    print {$fh} $info_of->{$species}{seq} . "\n";
                }
                close $fh;
            }
        }
        elsif ( $line =~ /^#/ ) {    # comments
            next;
        }
        elsif ( $line =~ /^s\s/ ) {    # s line, contain info and seq
            $content .= $line;
        }
        else {                         # a, i, e, q lines
                # just ommit it
                # see http://genome.ucsc.edu/FAQ/FAQformat.html#format5
        }
    }

    print ".fas file generated.\n\n";
};

# process each .fasta files
my $stopwatch = AlignDB::Stopwatch->new;

my @jobs = sort @files;
my $run  = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@jobs,
    code     => $worker,
);
$run->run;

$stopwatch->block_message( "All files have been processed.", "duration" );
exit;

__END__

=head1 NAME

    maf2fasta.pl - convert maf to fasta

=head1 SYNOPSIS
    perl maf2fasta.pl --in_dir G:/S288CvsRM11

    maf2fasta.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --in_dir            fasta files' location
        --out_dir           output location
        --length            length threshold
        --parallel          run in parallel mode

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
