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
use Bio::AlignIO;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);
use AlignDB::Run;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_dir = '.';    # Specify blocked fasta location here
my $file_vcf;        # Specify vcf location here
my $out_dir;         # Specify output dir here
my $nth = 1;         # which seq is target, start from 0

my $parallel = 1;    # run in parallel mode

my $verbose = 0;     # verbose mode

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'v|verbose'    => \$verbose,
    'i|in_dir=s'   => \$in_dir,
    'f|file_vcf=s' => \$file_vcf,
    'o|out_dir=s'  => \$out_dir,
    'n|nth=i'      => \$nth,
    'p|parallel=i' => \$parallel,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Preparing section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;

# make output dir
unless ($out_dir) {
    $out_dir = File::Spec->rel2abs($in_dir) . "_vcf";
}
if ( !-e $out_dir ) {
    mkdir $out_dir, 0777;
}
else {
    print "We are going to output blocked fasta.\n";
    die "$out_dir exists, you should remove it first to avoid errors.\n";
}

# Search for all files
my @files_fa = sort File::Find::Rule->file->name('*.fas')->in($in_dir);
printf "\n----Total .fa Files: %4s----\n\n", scalar @files_fa;

#----------------------------#
# read vcf
#----------------------------#
print "Parsing $file_vcf\n";
my ( $new_strain, $snvs, $indels, $others ) = parse_vcf($file_vcf);
print "Dumping other variations\n";
DumpFile( "others_vcf.yml", $others );
print "Splitting variations by chromosomes\n";
my $vars_of = {};
for my $var ( @{$snvs}, @{$indels} ) {
    my $chr = $var->[0];
    if ( exists $vars_of->{$chr} ) {
        push @{ $vars_of->{$chr} }, $var;
    }
    else {
        $vars_of->{$chr} = [$var];
    }
}
print "Sorting variations\n";
for my $chr ( sort keys %{$vars_of} ) {
    my @sorted = sort { $a->[1] <=> $b->[1] } @{ $vars_of->{$chr} };
    $vars_of->{$chr} = \@sorted;
}
print "Done.\n\n";

$stopwatch->block_message( ".vcf parsed.", "duration" );

#----------------------------------------------------------#
# run
#----------------------------------------------------------#
my $worker = sub {
    my $infile    = shift;
    my $stopwatch = AlignDB::Stopwatch->new;
    print "Process $infile\n";

    # don't use $/ = "\n\n", which cause bioperl panic
    open my $in_fh, "<", $infile;
    my $content = '';
READ: while ( my $line = <$in_fh> ) {
        if ( $line =~ /^\s+$/ and $content =~ /\S/ ) {
            my @lines = grep {/\S/} split /\n/, $content;
            $content = '';
            die "headers not equal to seqs\n" if @lines % 2;

            my ( $seq_of, $seq_names ) = ( {}, [] );
            while (@lines) {
                my $name = shift @lines;
                $name =~ s/^\>//;
                chomp $name;
                my $seq = shift @lines;
                chomp $seq;
                push @{$seq_names}, $name;
                $seq_of->{$name} = uc $seq;
            }

            # S288C.chrXIV(1):763984-765396|species=S288C
            my $target_str = $seq_names->[$nth];
            print "Processing $target_str\n";

            my $info_qr = qr{
                ([\w-]+)        # strain
                \.              # spacer
                (\w+)           # chr
                \((.{1,2})\)    # strand
                \:              # spacer
                (\d+)           # start
                \-              # spacer
                (\d+)           # end
            }x;

            if ( $target_str !~ $info_qr ) {
                print "Target info error: [$target_str]\n";
                next READ;
            }
            my $info = {
                target => $target_str,
                strain => $1,
                chr    => $2,
                strand => $3,
                start  => $4,
                end    => $5,
                seq    => $seq_of->{$target_str},
            };

            my $vars = $vars_of->{ $info->{chr} };

            # pass a smaller set of vars to insert_vars()
            # for more safe, expend the set a bit
            my $idx_start = find_pos( $vars, $info->{start} );
            for my $i ( 1 .. 10 ) {
                last if $idx_start < 1;
                $idx_start--;
            }
            my $idx_end = find_pos( $vars, $info->{end} );
            for my $i ( 1 .. 10 ) {
                last if $idx_end >= scalar @{$vars} - 1;
                $idx_end++;
            }
            my $smaller_vars = [ @{$vars}[ $idx_start .. $idx_end ] ];

            # insert vars
            my $new_seq = insert_vars( $info, $smaller_vars );

            my $gap_number = $new_seq =~ tr/-/-/;
            my $new_seq_length = length($new_seq) - $gap_number;
            my $new_name
                = "$new_strain.chrUn(+):1-$new_seq_length|species=$new_strain";
            push @{$seq_names}, $new_name;
            $seq_of->{$new_name} = $new_seq;

            my $outfile = basename($infile);
            $outfile = $out_dir . "/$outfile";

            open my $out_fh, '>>', $outfile;
            for my $name ( @{$seq_names} ) {
                print {$out_fh} ">", $name, "\n";
                print {$out_fh} $seq_of->{$name}, "\n";
            }
            print {$out_fh} "\n";
            close $out_fh;
        }
        else {
            $content .= $line;
        }
    }
    close $in_fh;
    print "Done.\n\n";
};

# process each .fas files
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@files_fa,
    code     => $worker,
);
$run->run;

$stopwatch->block_message( "All files have been processed.", "duration" );
exit;

#----------------------------------------------------------#
# subs
#----------------------------------------------------------#
sub parse_vcf {
    my $vcf_file = shift;

    my $strain;
    my $snvs   = [];
    my $indels = [];
    my $others = [];

    open my $fh, '<', $vcf_file;
    my $flag_strain;
    while ( my $line = <$fh> ) {
        chomp $line;
        if ( $line =~ /^#/ ) {
            if ( $line =~ /^#CHROM/ ) {
                my @fields = split /\t/, $line;
                if ( scalar @fields < 10 ) {
                    print "There are no strain name and no genotype column\n";
                    ($strain) = split /\W+/, basename($vcf_file);
                    $strain .= "_vcf";
                }
                else {
                    $strain = ( split /\.|\-/, ( split /\t/, $line )[9] )[0];
                    $flag_strain++;
                }
            }
            next;
        }

        my @fields   = split /\t/, $line;
        my $chr      = $fields[0];
        my $pos      = int( $fields[1] );
        my $ref      = $fields[3];
        my $alt      = $fields[4];
        my $genotype = $flag_strain ? $fields[9] : '1/1';

        my $array_ref = [ $chr, $pos, $ref, $alt, $genotype ];
        if ( $genotype !~ /1/ ) {
            next;
        }
        if ( length($ref) == 1 ) {
            if ( length($alt) == 1 ) {
                push @{$snvs}, $array_ref;
            }
            elsif ( $alt =~ /^[AGCT]+$/i ) {
                push @{$indels}, $array_ref;    # insertions
            }
            else {
                push @{$others}, $array_ref;
            }
        }
        elsif ( length($alt) == 1 ) {
            if ( $alt =~ /^[AGCT]+$/i ) {
                push @{$indels}, $array_ref;    # deletions
            }
            else {
                push @{$others}, $array_ref;
            }
        }
        else {
            push @{$others}, $array_ref;
        }
    }
    close $fh;

    return ( $strain, $snvs, $indels, $others );
}

sub insert_vars {
    my $info = shift;
    my $vars = shift;

    # seq info
    my $chr       = $info->{chr};
    my $chr_start = $info->{start};
    my $chr_end   = $info->{end};
    my $chr_set   = AlignDB::IntSpan->new("$chr_start-$chr_end");
    my $strand    = $info->{strand};
    my $seq       = $info->{seq};
    my $length    = length $seq;
    if ( $strand eq "-" or $strand eq "-1" ) {
        $seq = revcom($seq);
    }

    # to keep the original indel positions
    my $seq_set
        = AlignDB::IntSpan->new("1-$length")->diff( find_indel_set($seq) );
    print "seq set: ",    $seq_set->runlist, "\n" if $verbose;
    print "seq length: ", $seq_set->count,   "\n" if $verbose;

    # get variations belong to this region
    my @region_vars = grep { $chr_set->member( $_->[1] ) } @{$vars};

    # reverse order, so inserted indel will not affair posterior variations
    @region_vars = sort { $b->[1] <=> $a->[1] } @region_vars;
    print scalar @region_vars, " variations in this region\n";

    for my $var (@region_vars) {
        my $pos        = $var->[1];
        my $ref        = $var->[2];
        my $alt        = $var->[3];
        my $region_pos = $pos - $chr_start + 1;
        my $ref_length = length $ref;

        my @string_pos;
        for my $i ( 0 .. $ref_length - 1 ) {
            push @string_pos, $seq_set->at( $region_pos + $i );
        }

        my $chars;
        for my $i (@string_pos) {
            $chars .= substr $seq, $i - 1, 1;
        }

        if ( $chars ne $ref ) {
            print "ref chars are not identical\n";
            print Dump {
                region_pos => $region_pos,
                string_pos => \@string_pos,
                chars      => $chars,
                var        => $var,
            };
            print "\n";
        }
        else {
            print "replace [$ref] with [$alt] at $pos\n" if $verbose;

            my $rep_start  = $string_pos[0];
            my $rep_end    = $string_pos[-1];
            my $rep_length = $rep_end - $rep_start + 1;
            substr $seq, $rep_start - 1, $rep_length, $alt;
        }
    }

    if ( $strand eq "-" or $strand eq "-1" ) {
        $seq = revcom($seq);
    }
    return $seq;
}

# Return the index of the first element >= the supplied value.
sub find_pos {
    my $vars = shift;
    my $val  = shift;

    my $low  = 0;
    my $high = scalar( @{$vars} ) - 1;

    while ( $low < $high ) {
        my $mid = int( ( $low + $high ) / 2 );
        if ( $val < $vars->[$mid][1] ) {
            $high = $mid;
        }
        elsif ( $val > $vars->[$mid][1] ) {
            $low = $mid + 1;
        }
        else {
            return $mid;
        }
    }

    return $low;
}

__END__

=head1 NAME

    add_vcf2fa.pl - convert vcf to a new sequence in block fasta files

=head1 SYNOPSIS
    perl add_vcf2fa.pl  -i d:\data\alignment\yeast65\S288CvsXVIIIGE10m_mafft\ -f d:\data\yeast10000\ref\S288C.vcf

    add_vcf2fa.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --verbose           verbose mode
        -i, --in_dir        blocked fasta location here
        -f, --file_vcf      vcf file location
        -o, --out_dir       output location
        -n, --nth           which seq is target, start from 0, default is 1
        -p, --parallel      run in parallel mode

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut

XXX There're still bugs in insert_vars
