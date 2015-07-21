package MyUtil;
use strict;
use warnings;
use autodie;

use File::Spec;
use File::HomeDir;

use List::MoreUtils qw(firstidx uniq zip);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use WWW::Mechanize;
use HTML::TableExtract;
use Regexp::Common qw(balanced);

use Bio::DB::Taxonomy;
use Bio::TreeIO;

use base 'Exporter';
use vars qw(@ISA @EXPORT_OK %EXPORT_TAGS);
@ISA = qw(Exporter);

%EXPORT_TAGS = (
    all => [
        qw{
            replace_home wgs_worker find_ancestor find_group
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

sub wgs_worker {
    my $term = shift;

    my $mech = WWW::Mechanize->new;
    $mech->stack_depth(0);    # no history to save memory

    # local shadowsocks proxy
    if ($ENV{SSPROXY}) {
        $mech->proxy(['http','https'], 'socks://127.0.0.1:1080' )
    }

    my $url_part = "http://www.ncbi.nlm.nih.gov/Traces/wgs/";
    my $url      = $url_part . '?val=' . $term;
    print " " x 4 . $url . "\n";

    my $info = { prefix => $term, };
    $mech->get($url);

    {                         # extract from tables
        my $page    = $mech->content;
        my @tables  = qw{ meta-table structured-comments };
        my @columns = (
            '#_of_Contigs',    'Total_length',
            'Update_date',     'BioProject',
            'Keywords',        'Organism',
            'Assembly_Method', 'Assembly_Name',
            'Genome_Coverage', 'Sequencing_Technology',
            'Biosource',
        );

        for my $table (@tables) {
            print " " x 4 . "Extract from table ", $table, "\n";
            my $te
                = HTML::TableExtract->new( attribs => { class => $table, }, );
            $te->parse($page);

            my %srr_info;
            for my $ts ( $te->table_states ) {
                for my $row ( $ts->rows ) {
                    for my $cell (@$row) {
                        if ($cell) {
                            $cell =~ s/[,:]//g;
                            $cell =~ s/^\s+//g;
                            $cell =~ s/\s+$//g;
                            $cell =~ s/\s+/ /g;
                        }
                    }
                    next unless $row->[0];
                    $row->[0] =~ s/\s+/_/g;
                    next unless grep { $row->[0] eq $_ } @columns;

                    $row->[1] =~ s/\s+.\s+show.+lineage.+$//g;
                    if ( $row->[0] eq 'Biosource' ) {
                        my ($biosource_strain)
                            = grep {/strain = /} grep {defined} split /\//,
                            $row->[1];

                        #print $row->[1], "\n";
                        $biosource_strain =~ s/strain = //;

                        #print $biosource_strain, "\n";
                        $info->{ $row->[0] } = $biosource_strain;
                    }
                    else {
                        $info->{ $row->[0] } = $row->[1];
                    }
                }
            }
        }
    }

    {    # taxon id
        my @links = $mech->find_all_links( url_regex => => qr{wwwtax}, );
        if ( @links and $links[0]->url =~ /\?id=(\d+)/ ) {
            $info->{taxon_id} = $1;
        }
    }

    {    # pubmed id
        my @links = $mech->find_all_links( url_regex => => qr{\/pubmed\/}, );
        if ( @links and $links[0]->url =~ /\/pubmed\/(\d+)/ ) {
            $info->{pubmed} = $1;
        }
    }

    {    # downloads
        my @links = $mech->find_all_links(
            text_regex => qr{$term},
            url_regex  => => qr{download=},
        );
        $info->{download} = [ map { $url_part . $_->url } @links ];
    }

    return $info;
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
