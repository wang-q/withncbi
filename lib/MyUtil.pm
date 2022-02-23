package MyUtil;
use strict;
use warnings;
use autodie;

use File::HomeDir;
use Path::Tiny;

use List::MoreUtils qw(firstidx uniq zip);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use WWW::Mechanize;
use HTML::TableExtract;
use Regexp::Common qw(balanced);

use Bio::DB::Taxonomy;
use Bio::Taxon;
use Bio::TreeIO;

use base 'Exporter';
our (@ISA, @EXPORT_OK, %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            find_ancestor find_group abbr abbr_most
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

#@returns Bio::Taxon
sub find_ancestor {
    my Bio::Taxon $node = shift;
    my $rank = shift || 'species';

    return $node if $node->rank eq $rank;

RANK: while (1) {
        $node = $node->ancestor;
        last RANK unless defined $node;
        return $node if $node->rank eq $rank;
    }

    return;
}

### After remove 'no rank'
# chimpanzee
# superkingdom phylum subphylum class superorder order suborder infraorder parvorder superfamily family subfamily genus

# turkey
# superkingdom phylum subphylum class superorder order family subfamily genus

# rice
# superkingdom kingdom phylum class subclass order family subfamily tribe genus

# Salmonella
# superkingdom phylum class order family genus

# Methanococcus
# superkingdom phylum class order family genus

sub find_group {
    my $node = shift;

    my $tree = Bio::Tree::Tree->new( -node => $node );

    my @nodes = $tree->get_lineage_nodes($node);
    my @ranks = map { $_->rank } @nodes;
    my @names = map { $_->scientific_name } @nodes;

    my ( $group, $subgroup );
    my $unclassified;
    if ( grep {/unclassified|Candidatus|unspecified/i} @names ) {
        $unclassified = 1;
    }

    my $superkingdom = $names[ firstidx { $_ eq 'superkingdom' } @ranks ];

    if ($unclassified) {
        $group    = $superkingdom;
        $subgroup = 'unclassified';
    }
    else {
        my $kingdom_idx = firstidx { $_ eq 'kingdom' } @ranks;
        my $kingdom = $names[$kingdom_idx] if defined $kingdom_idx;
        my $phylum_idx = firstidx { $_ eq 'phylum' } @ranks;
        my $phylum = $names[$phylum_idx] if defined $phylum_idx;
        my $class_idx = firstidx { $_ eq 'class' } @ranks;
        my $class = $names[$class_idx] if defined $class_idx;

        if ( $superkingdom eq 'Bacteria' or $superkingdom eq 'Archaea' ) {
            if ( defined $phylum ) {
                $group = $phylum;
                if ( defined $class ) {
                    $subgroup = $class;
                }
                else {
                    $subgroup = 'Other';
                }
            }
            else {
                ( $group, $subgroup ) = ( $superkingdom, 'Other' );
            }
        }
        elsif ( $superkingdom eq 'Eukaryota' ) {

            # groups:
            # Animals, Fungi, Plants, Protists
            #
            # subgroup:
            # Mammals, Birds, Fishes, Flatworms, Insects, Amphibians, Reptiles,
            #   Roundworms
            # Ascomycetes, Basidiomycetes
            # Land Plants, Green Algae
            # Apicomplexans, Kinetoplasts

            if ( defined $kingdom ) {
                if ( $kingdom eq 'Metazoa' ) {
                    $group = 'Animals';

                    if ( defined $phylum ) {
                        if ( $phylum eq 'Platyhelminthes' ) {
                            $subgroup = 'Flatworms';
                        }
                        elsif ( $phylum eq 'Nematoda' ) {
                            $subgroup = 'Roundworms';
                        }
                        elsif ( $phylum eq 'Arthropoda' ) {
                            if ( grep { $_ eq 'Hexapoda' } @names ) {
                                $subgroup = 'Insects';
                            }
                        }
                        elsif ( $phylum eq 'Chordata' ) {
                            if ( grep {/^(Testudines|Lepidosauria|Crocodylia)$/}
                                @names )
                            {
                                $subgroup = 'Reptiles';
                            }
                            elsif (
                                grep {
                                    /^(Chondrichthyes|Dipnoi|Actinopterygii|Hyperotreti|Hyperoartia|Coelacanthimorpha)/
                                } @names
                                )
                            {
                                $subgroup = 'Fishes';
                            }
                            elsif ( defined $class ) {
                                if ( $class eq 'Mammalia' ) {
                                    $subgroup = 'Mammals';
                                }
                                elsif ( $class eq 'Aves' ) {
                                    $subgroup = 'Birds';
                                }
                                elsif ( $class eq 'Amphibia' ) {
                                    $subgroup = 'Amphibians';
                                }
                            }
                        }
                    }

                    if ( !defined $subgroup ) {
                        $subgroup = 'Other';
                    }
                }
                elsif ( $kingdom eq 'Fungi' ) {
                    $group = 'Fungi';
                    if ( defined $phylum ) {
                        if ( $phylum eq 'Ascomycota' ) {
                            $subgroup = 'Ascomycetes';
                        }
                        elsif ( $phylum eq 'Basidiomycota' ) {
                            $subgroup = 'Basidiomycetes';
                        }
                        else {
                            $subgroup = 'Other';
                        }
                    }
                }
                elsif ( $kingdom eq 'Viridiplantae' ) {
                    $group = 'Plants';
                    if ( grep { $_ eq 'Embryophyta' } @names ) {
                        $subgroup = 'Land Plants';
                    }
                    else {
                        $subgroup = 'Green Algae';
                    }
                }
                else {
                    $group = 'Protists';
                    if ( grep { $_ eq 'Apicomplexa' } @names ) {
                        $subgroup = 'Apicomplexans';
                    }
                    elsif ( grep { $_ eq 'Kinetoplastida' } @names ) {
                        $subgroup = 'Kinetoplasts';
                    }
                    else {
                        $subgroup = 'Other';
                    }
                }
            }
            else {
                ( $group, $subgroup ) = qw{Eukaryota Other};
            }
        }
        else {
            ( $group, $subgroup ) = qw{Other Other};
        }
    }

    return ( $superkingdom, $group, $subgroup );
}

# from core module Text::Abbrev
sub abbr {
    my $list = shift;
    my $min  = shift;
    return {} unless @{$list};

    $min = 1 unless $min;
    my $hashref = {};
    my %table;

WORD: for my $word ( @{$list} ) {
        for ( my $len = length($word) - 1; $len >= $min; --$len ) {
            my $abbrev = substr( $word, 0, $len );
            my $seen = ++$table{$abbrev};
            if ( $seen == 1 ) {    # We're the first word so far to have
                                   # this abbreviation.
                $hashref->{$abbrev} = $word;
            }
            elsif ( $seen == 2 ) {    # We're the second word to have this
                                      # abbreviation, so we can't use it.
                delete $hashref->{$abbrev};
            }
            else {                    # We're the third word to have this
                                      # abbreviation, so skip to the next word.
                next WORD;
            }
        }
    }

    # Non-abbreviations always get entered, even if they aren't unique
    for my $word ( @{$list} ) {
        $hashref->{$word} = $word;
    }
    return $hashref;
}

sub abbr_most {
    my $list  = shift;
    my $min   = shift;
    my $creat = shift; # don't abbreviate 1 letter. "I'd spell creat with an e."
    return {} unless @{$list};

    # don't abbreviate
    if ( defined $min and $min == 0 ) {
        my %abbr_of;
        for my $word ( @{$list} ) {
            $abbr_of{$word} = $word;
        }
        return \%abbr_of;
    }

    # do abbreviate
    else {
        my $hashref = $min ? abbr( $list, $min ) : abbr($list);
        my @ks = sort keys %{$hashref};
        my %abbr_of;
        for my $i ( reverse( 0 .. $#ks ) ) {
            if ( index( $ks[$i], $ks[ $i - 1 ] ) != 0 ) {
                $abbr_of{ $hashref->{ $ks[$i] } } = $ks[$i];
            }

            if ($creat) {
                for my $key ( keys %abbr_of ) {
                    if ( length($key) - length( $abbr_of{$key} ) == 1 ) {
                        $abbr_of{$key} = $key;
                    }
                }
            }
        }

        return \%abbr_of;
    }
}

1;
