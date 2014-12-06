package MyUtil;
use strict;
use warnings;
use autodie;

use File::Spec;
use File::HomeDir;

use List::MoreUtils qw(firstidx);

use Bio::DB::Taxonomy;
use Bio::TreeIO;

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            replace_home find_ancestor find_group
            },
    ],
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

sub replace_home {
    my $path = shift;

    if ( $path =~ /^\~\// ) {
        $path =~ s/^\~\///;
        $path = File::Spec->catdir( File::HomeDir->my_home, $path );
    }

    return $path;
}

sub find_ancestor {
    my $node = shift;
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

1;
