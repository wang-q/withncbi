# A simple class to get annotations for a certain chromosomal region (slice) using ensembl database
# and API.
#
# Qiang Wang

package MyEnsembl;
use Moose;
use Carp;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;

use AlignDB::IntSpan;

has 'server'     => ( is => 'ro', isa => 'Str' );
has 'db'         => ( is => 'ro', isa => 'Str' );
has 'port'       => ( is => 'ro', isa => 'Int', default => sub {3306} );
has 'user'       => ( is => 'ro', isa => 'Str' );
has 'passwd'     => ( is => 'ro', isa => 'Str' );
has 'db_adaptor' => ( is => 'ro', isa => 'Object' );
has 'slice_obj'  => ( is => 'ro', isa => 'Object' );
has 'slice'      => ( is => 'rw', isa => 'HashRef', default => sub { {} } );

sub BUILD {
    my $self = shift;

    # A hash of all attributes with default values
    $self->{slice} = {
        _chr            => 0,
        _start          => 0,
        _end            => 0,
        _slice_set      => '??',
        _intergenic_set => '??',
        _gene_set       => '??',
        _tc_set         => '??',
        _exon_set       => '??',
        _intron_set     => '??',
        _cds_set        => '??',
        _utr_set        => '??',
        _repeat_set     => '??',
        _non_repeat_set => '??',
        _te_set         => '??',
        _non_te_set     => '??',
        _ncnr_set       => '??',
    };

    # Connect to the ensembl database
    my $db_adaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -host   => $self->server,
        -dbname => $self->db,
        -port   => $self->port,
        -user   => $self->user,
        -pass   => $self->passwd,
    ) or Carp::confess "Cannot connect to EnsEMBL database\n";

    $self->{db_adaptor} = $db_adaptor;

    return;
}

sub new_slice_obj {
    my $self  = shift;
    my $chr   = shift;
    my $start = shift;
    my $end   = shift;

    # obtain a slice
    my Bio::EnsEMBL::DBSQL::DBAdaptor $db_adaptor       = $self->db_adaptor;
    my Bio::EnsEMBL::DBSQL::SliceAdaptor $slice_adaptor = $db_adaptor->get_SliceAdaptor;

    my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $end );

    return $slice;
}

sub set_slice {
    my $self  = shift;
    my $chr   = shift;
    my $start = shift;
    my $end   = shift;

    my $slice = $self->new_slice_obj( $chr, $start, $end );
    $self->{slice_obj} = $slice;

    # create feature sets
    my $slice_set = AlignDB::IntSpan->new("$start-$end");
    my $gene_set  = AlignDB::IntSpan->new;
    my $exon_set  = AlignDB::IntSpan->new;
    my $cds_set   = AlignDB::IntSpan->new;

    # Gene
    for my $gene ( @{ $slice->get_all_Genes } ) {
        $gene = $gene->transform('chromosome');
        $gene_set->merge( _ftr2runlist($gene) );
        for my $trans ( @{ $gene->get_all_Transcripts } ) {

            # Exon
            for my $exon ( @{ $trans->get_all_Exons } ) {
                $exon = $exon->transform('chromosome');
                $exon_set->merge( _ftr2runlist($exon) );
            }

            # CDS
            for my $cds ( @{ $trans->get_all_translateable_Exons } ) {
                $cds = $cds->transform('chromosome');
                $cds_set->merge( _ftr2runlist($cds) );
            }
        }
    }

    $gene_set = $gene_set->intersect($slice_set);
    $exon_set = $exon_set->intersect($slice_set);
    $cds_set  = $cds_set->intersect($slice_set);

    # others
    my $tc_runlist = '-';
    if ( $exon_set->is_not_empty ) {
        $tc_runlist = $exon_set->min . '-' . $exon_set->max;
    }
    my $tc_set         = AlignDB::IntSpan->new($tc_runlist);
    my $intergenic_set = $slice_set->diff($gene_set);
    my $intron_set     = $tc_set->diff($exon_set);
    my $utr_set        = $exon_set->diff($cds_set);

    # repeat set
    my $repeat_set = AlignDB::IntSpan->new;

    for my $repeat ( @{ $slice->get_all_RepeatFeatures } ) {
        $repeat = $repeat->transform('chromosome');
        $repeat_set->merge( _ftr2runlist($repeat) );
    }

    $repeat_set = $slice_set->intersect($repeat_set);

    my $non_repeat_set = $slice_set->diff($repeat_set);

    # transposon set
    my $te_set   = AlignDB::IntSpan->new;
    my @te_types = (
        "LTRs",
        "Type I Transposons\/LINE",
        "Type I Transposons\/SINE",
        "Type II Transposons",
        "TE",    # gramene type
    );

    for my $te_type (@te_types) {
        for my $te ( @{ $slice->get_all_RepeatFeatures( undef, $te_type ) } ) {
            $te = $te->transform('chromosome');
            $te_set->merge( _ftr2runlist($te) );
        }
    }
    $te_set = $slice_set->intersect($te_set);

    my $non_te_set = $slice_set->diff($te_set);

    my $ncnr_set = $slice_set->diff($cds_set)->intersect($non_repeat_set);

    # Store slice specified properties and feature sets
    $self->{slice}{_chr}            = $chr;
    $self->{slice}{_start}          = $start;
    $self->{slice}{_end}            = $end;
    $self->{slice}{_slice_set}      = $slice_set;
    $self->{slice}{_intergenic_set} = $intergenic_set;
    $self->{slice}{_gene_set}       = $gene_set;
    $self->{slice}{_exon_set}       = $exon_set;
    $self->{slice}{_tc_set}         = $tc_set;
    $self->{slice}{_intron_set}     = $intron_set;
    $self->{slice}{_cds_set}        = $cds_set;
    $self->{slice}{_utr_set}        = $utr_set;
    $self->{slice}{_repeat_set}     = $repeat_set;
    $self->{slice}{_non_repeat_set} = $non_repeat_set;
    $self->{slice}{_te_set}         = $te_set;
    $self->{slice}{_non_te_set}     = $non_te_set;
    $self->{slice}{_ncnr_set}       = $ncnr_set;

    return;
}

sub slice_stat {
    my $self = shift;

    my $slice_hash = $self->slice;

    my %stat_hash;

    for my $key ( keys %$slice_hash ) {
        my $new_key = $key;
        if ( $new_key =~ /_set/ ) {
            $new_key =~ s/_(.+)_set/$1/;
            $stat_hash{$new_key} = $slice_hash->{$key}->size;
        }
        else {
            $new_key =~ s/_//;
            $stat_hash{$new_key} = $slice_hash->{$key};
        }
    }

    return \%stat_hash;
}

sub locate_position {
    my ( $self, $pos_start, $pos_end ) = @_;

    # for single base, start equal to end
    $pos_end ||= $pos_start;

    if ( $pos_start > $pos_end ) {
        ( $pos_start, $pos_end ) = ( $pos_end, $pos_start );
    }

    my $pos_set = AlignDB::IntSpan->new("$pos_start-$pos_end");

    return $self->locate_set_position($pos_set);
}

sub locate_set_position {
    my $self = shift;
    my AlignDB::IntSpan $pos_set = shift;

    my $pos_location;      # location signal
    my $pos_in_repeats;    # repeats signal

    my $slice_hash = $self->slice;

    # Check $pos_start and pos_end are in the Slice
    unless ( $pos_set->subset( $slice_hash->{_slice_set} ) ) {
        Carp::carp "Range ", $pos_set->runlist, " is not in Slice ",
            $slice_hash->{_start}, "-", $slice_hash->{_end}, "!\n";
    }

    # one feature, three forms: coding; semi_coding; non_coding
    my @features = qw{
        _cds_set
    };

    for my $feature (@features) {
        if ( $pos_set->subset( $slice_hash->{$feature} ) ) {
            $pos_location = $feature;
        }
        else {
            my $intersect = $pos_set->intersect( $slice_hash->{$feature} );
            my $n         = $intersect->size;
            if ( $n >= 1 ) {
                $pos_location = "semi_" . $feature;
            }
            else {
                $pos_location = "non_" . $feature;
            }
        }
    }
    $pos_location =~ s/_cds_set/coding/;

    # one feature, three forms: repeat; semi_repeat; non_repeat
    my @repeats = qw{
        _repeat_set
    };

    for my $repeat (@repeats) {
        if ( $pos_set->subset( $slice_hash->{$repeat} ) ) {
            $pos_in_repeats = $repeat;
        }
        else {
            my $intersect = $pos_set->intersect( $slice_hash->{$repeat} );
            my $n         = $intersect->size;
            if ( $n >= 1 ) {
                $pos_in_repeats = "semi_" . $repeat;
            }
            else {
                $pos_in_repeats = "non_" . $repeat;
            }
        }
    }
    $pos_in_repeats =~ s/_(.+)_set/$1/;

    return ( $pos_location, $pos_in_repeats );
}

sub feature_portion {

    # start, end should be "$start-$end" format
    my ( $self, $feature, $pos_set ) = @_;

    my $slice_hash = $self->slice;

    my $set;
    if ( ref $pos_set eq "AlignDB::IntSpan" ) {
        $set = $pos_set;
    }
    else {
        $set = AlignDB::IntSpan->new($pos_set);
    }

    my $pos_n = $set->size;
    return if $pos_n <= 0;
    my $intersect       = $set->intersect( $slice_hash->{$feature} );
    my $n               = $intersect->size;
    my $feature_portion = $n / $pos_n;

    return $feature_portion;

}

# return an AlignDB::IntSpan object of a feature
sub feature_set_obj {
    my $self    = shift;
    my $feature = shift;

    my $slice_hash = $self->slice;

    return $slice_hash->{$feature};
}

sub _ftr2runlist {
    my $f     = shift;
    my $start = $f->start;
    my $end   = $f->end;
    return "$start-$end";
}

1;
