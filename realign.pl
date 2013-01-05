#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Basename;
use Parallel::ForkManager;

use AlignDB::Stopwatch;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $dir_align        = 'F:/S288CvsRM11'; # Specify location here
my $dir_out          = undef;            # Specify output dir here
my $length_threshold = 5000;             # Get the threshold of alignment length
my $aln_prog         = 'clustalw';       # Default alignment program

my $quick_mode   = undef;    # quick mode
my $expand_indel = 1000;     # in quick mode, expand indel regoin
my $join_indel   = 500;      # in quick mode, join adjacent indel regions

# run in parallel mode
my $parallel = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'         => \$help,
    'man'            => \$man,
    'da|dir_align=s' => \$dir_align,
    'do|dir_out=s'   => \$dir_out,
    'lt|length=i'    => \$length_threshold,
    'msa=s'          => \$aln_prog,
    'quick'          => \$quick_mode,
    'expand=i'       => \$expand_indel,
    'join=i'         => \$join_indel,
    'parallel=i'     => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# make output dir
#----------------------------------------------------------#
unless ($dir_out) {
    $dir_out = $dir_align . "_$aln_prog";
    $dir_out = $dir_out . "_quick" if $quick_mode;
}
unless ( -e $dir_out ) {
    mkdir $dir_out, 0777
        or die "Cannot create \"$dir_out\" directory: $!";
}

#----------------------------------------------------------#
# Search for all files and push their paths to @dir_fiels
#----------------------------------------------------------#
my @axt_files = File::Find::Rule->file()->name('*.axt')->in($dir_align);
printf "\n----Total .AXT Files: %4s----\n\n", scalar @axt_files;

if ( scalar @axt_files < $parallel ) {
    $parallel = scalar @axt_files;
}

#----------------------------------------------------------#
# realign
#----------------------------------------------------------#
my $process_axt = sub {
    my $infile     = shift;
    my $start_time = time;

    my $stopwatch = AlignDB::Stopwatch->new();
    $stopwatch->block_message("Process $infile...");

    my $outfile = basename($infile);
    $outfile = $dir_out . "/$outfile";

    open my $axt_fh, '<', $infile
        or die("Cannot open IN file $infile");
    open my $out_fh, '>', $outfile
        or die("Cannot open OUT file $outfile");

    while (1) {
        my $summary_line = <$axt_fh>;
        unless ($summary_line) {
            last;
        }
        if ( $summary_line =~ /^#/ ) {
            next;
        }
        chomp $summary_line;
        chomp( my $first_line = <$axt_fh> );
        $first_line = uc $first_line;
        chomp( my $second_line = <$axt_fh> );
        $second_line = uc $second_line;
        my $dummy = <$axt_fh>;

        my $align_length = length $first_line;
        unless ( $align_length > $length_threshold ) {
            next;
        }

        if ($quick_mode) {    # only realign regions containing indels

            # expand indel set
            my $first_indel_set = &find_indel_set( $first_line, $expand_indel );
            my $second_indel_set
                = &find_indel_set( $second_line, $expand_indel );

            my $realign_region = $first_indel_set->union($second_indel_set);
            $realign_region = $realign_region->join_span($join_indel);

            # realign all segments in realign_region
            # realign backward, so coordinate changs will not affect anything
            my @realign_region_spans = $realign_region->spans();
            foreach ( reverse @realign_region_spans ) {
                my $seg_start = $_->[0];
                my $seg_end   = $_->[1];
                my @segments;
                foreach ( $first_line, $second_line ) {
                    my $seg = substr( $_, $seg_start - 1,
                        $seg_end - $seg_start + 1 );
                    push @segments, $seg;
                }
                my $realign_segments = &multi_align( \@segments, $aln_prog );
                foreach ( $first_line, $second_line ) {
                    my $seg = shift @$realign_segments;
                    substr(
                        $_,
                        $seg_start - 1,
                        $seg_end - $seg_start + 1, $seg
                    );
                }
            }
        }
        else {
            my $realigned
                = &multi_align( [ $first_line, $second_line ], $aln_prog );

            $first_line  = $realigned->[0];
            $second_line = $realigned->[1];
        }

        print {$out_fh} $summary_line, "\n";
        print {$out_fh} $first_line,   "\n";
        print {$out_fh} $second_line,  "\n";
        print {$out_fh} "\n";
    }

    close $out_fh;
    close $axt_fh;

    $stopwatch->block_message( "$infile has been processed.", "duration" );
};

# process each .axt files
my $stopwatch = AlignDB::Stopwatch->new();
if ( $parallel > 1 ) {
    my $pm = new Parallel::ForkManager($parallel);

    foreach my $infile ( sort @axt_files ) {
        $pm->start and next;
        &{$process_axt}($infile);
        $pm->finish;
    }
    $pm->wait_all_children;
}
else {
    foreach my $infile ( sort @axt_files ) {
        &{$process_axt}($infile);
    }
}

$stopwatch->block_message( "All files have been processed.", "duration" );
exit;

__END__

=head1 NAME

    realign.pl - realign axt file

=head1 SYNOPSIS
    perl realign.pl -a=G:/S288CvsRM11 --msa=muscle --quick --expand=1000

    realign.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --dir_align         .axt files' directory
        --dir_out           output directory
        --length            length threshold
        --msa               alignment program
        --quick             use quick mode   
        --expand            in quick mode, expand indel region
        --join              in quick mode, join adjacent indel regions
        --parallel          run in parallel mode

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
