package MyKaKs;
use Moose;
use Carp;

use Bio::Seq;

# Ka/Ks estimators
use Bio::Tools::Run::Phylo::PAML::Codeml;

# Multiple Sequence Alignment programs
use Bio::Tools::Run::Alignment::Clustalw;

# for projecting alignments from protein to R/DNA space
use Bio::Align::Utilities qw(aa_to_dna_aln);

use App::Fasops::Common;

has 'aln_factory'  => ( is => 'ro', isa => 'Object', );
has 'kaks_factory' => ( is => 'ro', isa => 'Object', );

# fasta filename
has 'fasta' => ( is => 'rw', isa => 'Str', );

# seqs hashref
has 'seq_of' => ( is => 'rw', isa => 'HashRef', );

# results
has 'results' => ( is => 'ro', isa => 'ArrayRef', );

sub BUILD {
    my $self = shift;

    if ( $self->fasta ) {
        my $seq_of =  App::Fasops::Common::read_fasta( $self->fasta );
        $self->seq_of($seq_of);
    }

    my $aln_factory
        = Bio::Tools::Run::Alignment::Clustalw->new( -verbose => -1 );
    unless ( $aln_factory->executable ) {
        Carp::confess "Could not find the executable for clustalw\n";
    }
    $self->{aln_factory} = $aln_factory;

    my $kaks_factory = Bio::Tools::Run::Phylo::PAML::Codeml->new(
        -params => {
            'runmode' => -2,
            'seqtype' => 1,
        },
        -verbose => -1,
    );
    unless ( $kaks_factory->executable ) {
        Carp::confess "Could not find the executable for codeml\n";
    }
    $self->{kaks_factory} = $kaks_factory;

    return;
}

sub run {
    my $self = shift;

    my $dna_of = {};
    my @prots;
    for my $header ( keys %{ $self->seq_of } ) {
        my $seq = Bio::Seq->new(
            -display_id => $header,
            -seq        => $self->seq_of->{$header},
        );
        $dna_of->{$header} = $seq;

        my Bio::Seq $protein = $seq->translate;

        # check seqs
        my $pseq = $protein->seq;
        if ( $pseq =~ /\*/ and $pseq !~ /\*$/ ) {
            Carp::confess "provided a cDNA ["
                . $seq->display_id
                . "] sequence with stop codon(s), PAML will choke!\n";
        }

        push @prots, $protein;
    }

    my $aa_aln  = $self->aln_factory->align( \@prots );
    my $dna_aln = aa_to_dna_aln( $aa_aln, $dna_of );
    my @each    = $dna_aln->each_seq;
    $self->kaks_factory->alignment($dna_aln);

    my ( $rc, $parser ) = $self->kaks_factory->run;
    if ( $rc <= 0 ) {
        Carp::confess $self->kaks_factory->error_string, "\n";
    }
    my $result   = $parser->next_result;
    my $MLmatrix = $result->get_MLmatrix;

    my @otus = $result->get_seqs;
    my @pos  = map {
        my $c = 1;
        for my $s (@each) {
            last if ( $s->display_id eq $_->display_id );
            $c++;
        }
        $c;
    } @otus;

    my @results;
    for ( my $i = 0; $i < ( scalar @otus - 1 ); $i++ ) {
        for ( my $j = $i + 1; $j < ( scalar @otus ); $j++ ) {
            my $sub_aa_aln = $aa_aln->select_noncont( $pos[$i], $pos[$j] );
            my $sub_dna_aln = $dna_aln->select_noncont( $pos[$i], $pos[$j] );
            push @results,
                [
                $otus[$i]->display_id,
                $otus[$j]->display_id,
                $MLmatrix->[$i]->[$j]->{'dN'},
                $MLmatrix->[$i]->[$j]->{'dS'},
                $MLmatrix->[$i]->[$j]->{'omega'},
                $sub_aa_aln->percentage_identity,
                $sub_dna_aln->percentage_identity,
                ];
        }
    }
    $self->{results} = \@results;

    return;
}

1;
