#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

my $file_target;
my $file_merge;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?' => \$help,
    'man'    => \$man,
    't|ft=s' => \$file_target,
    'm|fm=s' => \$file_merge,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

for ( $file_target, $file_merge ) {
    die "Can't find file [$_]\n" unless -e $_;
}

#----------------------------------------------------------#
# apply
#----------------------------------------------------------#
my ( $count_t, $count_m ) = ( 0, 0 );

# read in merge
my $id_m   = [];
my $line_m = [];
{
    open my $fh, '<', $file_merge;
    while (<$fh>) {
        chomp;
        next unless $_;
        push @{$id_m}, ( split /,/ )[0];
        push @{$line_m}, $_;
        $count_m++;
    }
    close $fh;
}

# read in target
my $line_t = [];
{
    open my $fh, '<', $file_target;
    while (<$fh>) {
        chomp;
        next unless $_;
        my $id = ( split /,/ )[0];
        if ( grep { $_ eq $id } @{$id_m} ) {
            $count_t++;
            next;
        }
        else {
            push @{$line_t}, $_;
        }
    }
    close $fh;
}

# write target
{
    push @{$line_t}, @{$line_m};
    open my $fh, '>', $file_target;
    for ( @{$line_t} ) {
        print {$fh} $_, "\n";
    }
    close $fh;
}

print "Remove [$count_t] lines from $file_target\n";
print "Add [$count_m] lines from $file_merge\n";

exit;

__END__


=head1 NAME

    merge_csv.pl - Merge two csv files based on the first filed

=head1 SYNOPSIS

    merge_csv.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -t, --ft            target file (output)
        -m, --fm            merge file

    perl merge_csv.pl -t ../init/taxon.csv -m d:\data\alignment\yeast_combine\taxon.csv
    perl merge_csv.pl -t ../init/chr_length.csv -m d:\data\alignment\yeast_combine\chr_length.csv

=cut
