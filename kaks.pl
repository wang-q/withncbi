#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Spec;
use Math::Combinatorics;

use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Phylo::PAML::Yn00;
use Bio::Align::PairwiseStatistics;
use Bio::Align::Utilities qw(aa_to_dna_aln);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $input = 'cdna.fasta';
my $output;

my $localtmp = '';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'input=s'   => \$input,
    'output=s'  => \$output,
    'localtmp!' => \$localtmp,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

$output ||= $input . ".csv";

#----------------------------------------------------------#
# overide bioperl Clustalw wrapper tempdir, codeml tempdir is OK
#----------------------------------------------------------#
my $tmpdir;
{
    no warnings 'all';
    my $tempdir = sub {
        if ( defined $tmpdir ) {
            mkdir $tmpdir if !-e $tmpdir;
            return $tmpdir;
        }
        $tmpdir = 'tmp' . int( rand() * 10_000_000 );
        if ($localtmp) {
            $tmpdir = File::Spec->catdir( File::Spec->curdir(), $tmpdir );
        }
        else {
            $tmpdir = File::Spec->catdir( File::Spec->tmpdir(), $tmpdir );
        }
        mkdir $tmpdir if !-e $tmpdir;
        return $tmpdir;
    };
    *Bio::Tools::Run::Alignment::Clustalw::tempdir = $tempdir;
    use warnings;
}

#----------------------------------------------------------#
# Init align objects
#----------------------------------------------------------#
my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new;

my $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new(
    -params => {
        'runmode' => -2,
        'seqtype' => 1,
    }
);

unless ( $aln_factory->executable ) {
    die "Could not find the executable for clustalw\n";
}
unless ( $kaks_factory->executable ) {
    die "Could not find the executable for codeml\n";
}

#----------------------------------------------------------#
# read out seqs
#----------------------------------------------------------#
my ( $headers, $sequence_of ) = open_fasta($input);
if ( @$headers > keys %$sequence_of ) {
    die "There are duplicated sequence names in input file\n";
}

my $dna_of = {};
my $pro_of = {};
for my $header (@$headers) {
    my $seq = Bio::Seq->new(
        -display_id => $header,
        -seq        => $sequence_of->{$header},
    );
    $dna_of->{$header} = $seq;

    $pro_of->{$header} = $seq->translate;

    # check seqs
    my $pseq = $pro_of->{$header}->seq;
    if ( $pseq =~ /\*/ and $pseq !~ /\*$/ ) {
        die "provided a cDNA ["
            . $seq->display_id
            . "] sequence with stop codon(s), PAML will choke!\n";
    }
}

open my $out_fh, ">", $output or die "cannot open output file: $!\n";
print {$out_fh} join(
    ",",
    qw(Seq1 Seq1_Length Seq2 Seq2_Length
        Ka Ks Ka/Ks PI Prot_PercentID cDNA_PercentID)
    ),
    "\n";

my $combi = Math::Combinatorics->new(
    count => 2,
    data  => $headers,
);

print "============================\n";
print "Start at ", scalar localtime, "\n";
my $start_time = time;

my $comp_count = 0;
while ( my @combi_cdna = $combi->next_combination ) {

    my %seqs  = map { $_ => $dna_of->{$_} } @combi_cdna;
    my @dnas  = @{$dna_of}{@combi_cdna};
    my @prots = @{$pro_of}{@combi_cdna};

    my $raw_dna_aln = $aln_factory->align( \@dnas );
    my $stat        = Bio::Align::PairwiseStatistics->new;

    my $bases = $stat->number_of_comparable_bases($raw_dna_aln);
    my $nd    = $stat->number_of_differences($raw_dna_aln);
    my $ng    = $stat->number_of_gaps($raw_dna_aln);
    my $pi    = $nd / $bases;

    my $aa_aln  = $aln_factory->align( \@prots );
    my $dna_aln = aa_to_dna_aln( $aa_aln, \%seqs );
    my @each    = $dna_aln->each_seq;

    $kaks_factory->alignment($dna_aln);

    my ( $rc, $parser ) = $kaks_factory->run;
    if ( $rc <= 0 ) {
        die $kaks_factory->error_string, "\n";
    }
    my $result   = $parser->next_result;
    my $MLmatrix = $result->get_MLmatrix;

    my @otus = $result->get_seqs;

    print {$out_fh} join( ",",
        $dnas[0]->display_id,
        $dnas[0]->length,
        $dnas[1]->display_id,
        $dnas[1]->length,
        $MLmatrix->[0]->[1]->{'dN'},
        $MLmatrix->[0]->[1]->{'dS'},
        $MLmatrix->[0]->[1]->{'omega'},
        $pi,
        sprintf( "%.2f", $aa_aln->percentage_identity ),
        sprintf( "%.2f", $dna_aln->percentage_identity ),
        ),
        "\n";

    $comp_count++;
}

print "\n";
print "Total cDNA ", scalar @$headers, "\n";
print "Total combinations ", combi_k_n( 2, scalar @$headers ), "\n";
print "Comp\'ed combinations ", $comp_count, "\n";
print "\n";

print "End   at ", scalar localtime, "\n";
my $end_time  = time;
my $cost_time = $end_time - $start_time;
print "Cost $cost_time seconds.\n";
print "============================\n\n";

#----------------------------------------------------------#
# subroutines
#----------------------------------------------------------#
sub combi_k_n {
    my ( $k, $n ) = @_;
    my $com = factorial($n) / factorial($k) / factorial( $n - $k );
    return $com;
}

sub open_fasta {
    my $filename = shift;

    open my $in_fh, "<", $filename or die "cannot open input file: $!\n";
    my @contents = <$in_fh>;
    close $in_fh;

    my @headers;
    my %sequence_of;
    foreach my $current_line (@contents) {
        if ( $current_line =~ /^\>([\w-])+/ ) {
            chomp $current_line;
            $current_line =~ s/\>//;
            push @headers, $current_line;
            $sequence_of{$current_line} = '';
        }
        elsif ( $current_line =~ /^[\w-]+/ ) {
            chomp $current_line;
            my $header = $headers[-1];
            $sequence_of{$header} .= $current_line;
        }
        else {    # Blank line, do nothing.
        }
    }
    @headers = sort @headers;

    return ( \@headers, \%sequence_of );
}

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
