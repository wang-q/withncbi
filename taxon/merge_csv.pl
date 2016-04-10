#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use List::MoreUtils qw(firstidx);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $file_target;
my $file_merge;

my @fields;
my @fields2;

my $concat;
my $stdout;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    't|ft=s'     => \$file_target,
    'm|fm=s'     => \$file_merge,
    'f|fields=s' => \@fields,
    'f2|fields2=s' => \@fields2,
    'concat'     => \$concat,
    'stdout'     => \$stdout,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

for ( $file_target, $file_merge ) {
    die "Can't find file [$_]\n" unless -e $_;
}

if ( !scalar @fields ) {
    @fields = (0);
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
        my $id;
        if ($concat and scalar @fields2) {
            $id = join( "_", ( split /,/ )[@fields2] );
        }
        else {
            $id = join( "_", ( split /,/ )[@fields] );
        }
        push @{$id_m}, $id;
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
        my $id = join( "_", ( split /,/ )[@fields] );
        if ( grep { $_ eq $id } @{$id_m} ) {
            $count_t++;
            my $idx = firstidx { $_ eq $id } @{$id_m};
            splice @{$id_m}, $idx, 1;
            my ($line) = splice @{$line_m}, $idx, 1;

            my $newline;
            if ($concat) {
                $newline = "$_,$line";
            }
            else {
                $newline = $line;
            }
            push @{$line_t}, $newline;
        }
        else {
            push @{$line_t}, $_;
        }
    }
    close $fh;
}

#----------------------------#
# write outputs
#----------------------------#
if ( !$concat ) {
    push @{$line_t}, @{$line_m};
}

if ($stdout) {
    for ( @{$line_t} ) {
        print $_, "\n";
    }
}
else {
    open my $fh, '>', $file_target;
    for ( @{$line_t} ) {
        print {$fh} $_, "\n";
    }
    close $fh;
}

printf STDERR "%s [%s] lines from %s\n", $concat ? "Concat" : "Remove",
    $count_t,
    $file_target;
printf STDERR "Add [%s] lines from %s\n", $count_m, $file_merge;

exit;

__END__


=head1 NAME

merge_csv.pl - Merge two csv files based on @fields

=head1 SYNOPSIS

    merge_csv.pl [options]
      Options:
        --help                  brief help message
        --man                   full documentation
        -t, --ft STR            target file (output)
        -m, --fm SRT            merge file
        -f, --fields @INT       fields as identifies
        --concat                do concat other than merge
        -f2, --fieldss @INT     fields in concat file
        --stdout                write outputs to stdout instead of the target file

    perl merge_csv.pl -t ../init/taxon.csv -m d:\data\alignment\yeast_combine\taxon.csv 
    perl merge_csv.pl -t ../init/chr_length.csv -m d:\data\alignment\yeast_combine\chr_length.csv

=cut
