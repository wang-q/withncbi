#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use YAML::Syck;

use Path::Tiny;

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
        --concat    -c          do concat other than merge. Keep first ID fields

    cat 1.csv 2.csv | perl merge_csv.pl -f 0 -f 1

=cut

GetOptions(
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'outfile|o=s' => \( my $outfile = 'stdout' ),
    'fields|f=s'  => \my @fields,
    'concat|c'    => \my $concat,
) or Getopt::Long::HelpMessage(1);

if ( !scalar @fields ) {
    @fields = (0);
}
@fields = sort @fields;    # make array splicing happier

#----------------------------------------------------------#
# apply
#----------------------------------------------------------#

my $index_of = {};    # index of ids in @lines
my @lines;
my ( $count_all, $index ) = ( 0, 0 );

while ( my $line = <> ) {
    chomp $line;
    next unless $line;
    $count_all++;
    my $id = join( "_", ( split ",", $line )[@fields] );
    if ( exists $index_of->{$id} ) {
        if ($concat) {
            my $ori_index = $index_of->{$id};
            my $ori_line  = $lines[$ori_index];

            my @fs = split ",", $line;
            for my $f_idx ( reverse @fields ) {
                splice @fs, $f_idx, 1;
            }
            $lines[$ori_index] = join ",", $ori_line, @fs;
        }
    }
    else {
        $index_of->{$id} = $index;
        push @lines, $line;
        $index++;
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
        warn YAML::Syck::Dump { fields => \%seen, };
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

printf STDERR "Total lines [%d]; Result lines [%d].\n", $count_all, scalar @lines;

exit;

__END__
