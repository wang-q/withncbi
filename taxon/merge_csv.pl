#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use List::MoreUtils qw(firstidx);
use Path::Tiny;
use Set::Scalar;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

merge_csv.pl - Merge csv files based on @fields

=head1 SYNOPSIS

    cat input.csv another.csv | perl merge_csv.pl [options]
      Options:
        --help      -?          brief help message
        --outfile   -o  STR     output filename. Default is [stdout] for screen
        --fields    -f  @INT    fields as identifies. Default is [0], first column
        --concat                do concat other than merge. Keep first ID fields

    cat 1.csv 2.csv | perl merge_csv.pl -f 0 -f 1

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'outfile|o=s' => \( my $outfile = 'stdout' ),
    'fields|f=s'  => \my @fields,
    'concat'      => \my $concat,
) or HelpMessage(1);

if ( !scalar @fields ) {
    @fields = (0);
}
@fields = sort @fields;    # make array splicing happier

#----------------------------------------------------------#
# apply
#----------------------------------------------------------#
my ( $count_all, $count_op ) = ( 0, 0 );

my $id_set = Set::Scalar->new;
my ( @ids, @lines );

while ( my $line = <> ) {
    chomp $line;
    next unless $line;
    $count_all++;
    my $id = join( "_", ( split ",", $line )[@fields] );
    if ( $id_set->has($id) ) {
        $count_op++;
        if ($concat) {
            my $index = firstidx { $_ eq $id } @ids;
            my $ori_line = $lines[$index];

            my @fs = split ",", $line;
            for my $f_idx ( reverse @fields ) {
                splice @fs, $f_idx, 1;
            }
            $lines[$index] = join ",", $ori_line, @fs;
        }
    }
    else {
        $id_set->insert($id);
        push @ids,   $id;
        push @lines, $line;
    }
}

#----------------------------#
# check
#----------------------------#
{
    my %seen;
    for (@lines) {
        my $number = scalar split(",");
        $seen{$number}++;
    }
    if ( keys(%seen) > 1 ) {
        warn "*** Fields not identical, be careful.\n";
        warn Dump { fields => \%seen, };
    }
}

#----------------------------#
# write outputs
#----------------------------#
my $out_fh;
if ( lc($outfile) eq "stdout" ) {
    $out_fh = *STDOUT;
}
else {
    open $out_fh, ">", $outfile;
}

for (@lines) {
    print {$out_fh} $_ . "\n";
}
close $out_fh;

printf STDERR "Total lines [%d]; Result lines [%d]. ", $count_all, scalar @lines;
printf STDERR "%s [%d] lines.\n", $concat ? "Concat" : "Merge", $count_op;

exit;

__END__
