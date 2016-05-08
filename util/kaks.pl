#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::KaKs;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $input = 'cdna.fasta';
my $output;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'   => \$help,
    'man'      => \$man,
    'input=s'  => \$input,
    'output=s' => \$output,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$output ||= $input . ".csv";

#----------------------------------------------------------#
# run
#----------------------------------------------------------#
my $kaks = AlignDB::KaKs->new(fasta => $input);
$kaks->run;

open my $out_fh, ">", $output or die "cannot open output file: $!\n";
print {$out_fh}
    join( ",", qw{ Seq1 Seq2 Ka Ks Ka/Ks Prot_PercentID cDNA_PercentID } ),
    "\n";

for my $result (@{$kaks->results}) {
    print {$out_fh} join( ",", @$result), "\n";
}

close $out_fh;

__END__

=head1 NAME

    kaks.pl - Calculate Ka/Ks for cDNA

=head1 SYNOPSIS

    kaks.pl [options]
        Options:
            --help          brief help message
            --man           full documentation
            --input         input cdna filename
            --output        output filename
            --localtmp      use cwd as temporary dir

=cut
