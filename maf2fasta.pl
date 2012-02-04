#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Spec;
use File::Find::Rule;
use Bio::AlignIO;

use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);
use AlignDB::Run;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_dir = '.';    # Specify location here
my $out_dir;         # Specify output dir here
my $length_threshold = 1000;    # Set the threshold of alignment length
my $target_id        = 9606;    # taxon id of human
my $has_outgroup;    # Treate last sequence as outgroup and move it to top
my $subset;          # get sequences of listed names, seperated by comma

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'i|in_dir=s'   => \$in_dir,
    'o|out_dir=s'  => \$out_dir,
    'length=i'     => \$length_threshold,
    'id=i'         => \$target_id,
    'has_outgroup' => \$has_outgroup,
    'subset=s'     => \$subset,
    'parallel=i'   => \$parallel,
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
    $out_dir = $out_dir . "_$subset" if $subset;
    $out_dir = $out_dir . "_ref_first" if $has_outgroup;
}
unless ( -e $out_dir ) {
    mkdir $out_dir, 0777
        or die "Cannot create [$out_dir] directory: $!";
}

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my @files = File::Find::Rule->file()->name('*.maf')->in($in_dir);
printf "\n----Total .maf Files: %4s----\n\n", scalar @files;

#----------------------------------------------------------#
# run
#----------------------------------------------------------#
my $worker = sub {
    my $infile = shift;

    my $in = Bio::AlignIO->new(
        -file   => $infile,
        -format => 'maf'
    );

ALN: while ( my $aln = $in->next_aln ) {

        my $align_length  = $aln->length;
        my $num_sequences = $aln->num_sequences;

        next if $align_length < $length_threshold;

        my @names;
        my $seq_of = {};
        foreach my $seq ( $aln->each_seq ) {
            my $seq_id = $seq->display_id;
            my ($species) = split /\./, $seq_id;
            push @names, $species;
            $seq_of->{$species} = $seq;
        }

        if ($subset) {
            my $name_str = join " ", @names;
            my @subsets = split ",", $subset;
            next ALN if $num_sequences < @subsets;
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
        my ( undef, $chr ) = split /\./, $seq_of->{$target}->display_id;
        my $start = $seq_of->{$target}->start;
        my $end   = $seq_of->{$target}->end;

        open my $fh, '>',
            $out_dir . "/id3702" . "_$chr" . "_$start" . "_$end.fas";
        for (@names) {
            print {$fh} ">$_\n";
            print {$fh} $seq_of->{$_}->seq, "\n";
        }
        close $fh;
    }
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
