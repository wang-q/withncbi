#!/usr/bin/perl
use strict;
use warnings;

use Template;
use File::Basename;
use File::Find::Rule;
use File::Remove qw(remove);
use File::Spec;
use String::Compare;
use YAML qw(Dump Load DumpFile LoadFile);

#----------------------------------------------------------#
# RepeatMasker has been done
#----------------------------------------------------------#
my $store_dir = shift
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast_combine" );
my $parallel = 12;

{    # on linux
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast_combine" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    my @data = (

        # the best
        { taxon => 226125, name => 'Spar',   coverage => '7x', },
        { taxon => 285006, name => 'RM11',   coverage => '10x', },
        { taxon => 307796, name => 'YJM789', coverage => '10x', },

        # high coverage
        {   taxon    => 574961,
            name     => 'JAY291',
            coverage => '12x 454; 58x solexa se; 95x solexa pe',
        },
        {   taxon    => 538975,
            name     => 'Sigma1278b',
            coverage => '45x sanger/solexa',
        },
        {   taxon    => 643680,
            name     => 'EC1118',
            coverage => '17.6x 454; 4.1+1.9x sanger',
        },

        {   taxon    => 721032,
            name     => 'Kyokai_no__7',
            coverage => '9.1x sanger',
        },

        # wustl 11 yeast strains
        { taxon => 929585, name => 'T7', coverage => '25.4x 454/sanger', },

        # sgrp2 > 30x
        { taxon => 901501, name => "YPS128",    coverage => "56x illumina", },
        { taxon => 901503, name => "SK1",       coverage => "40x illumina", },
        { taxon => 901504, name => "L1528",     coverage => "34x illumina", },
        { taxon => 901506, name => "UWOPS83",   coverage => "638x illumina", },
        { taxon => 901507, name => "UWOPS87",   coverage => "834x illumina", },
        { taxon => 901509, name => "W303",      coverage => "33x illumina", },
        { taxon => 901511, name => "DBVPG6765", coverage => "37x illumina", },
        { taxon => 901512, name => "DBVPG6044", coverage => "32x illumina", },
    );

    for my $item ( sort @data ) {
        my $name = $item->{name};

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $name );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    print Dump \@data;

    my $tt = Template->new;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Saccharomyces,cerevisiae,[% item.name %],[% item.name %],
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "taxon.csv" )
    ) or die Template->error;

    # chr_length_chrUn.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]/WGS/sgrp2
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "chr_length_chrUn.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]/WGS/sgrp2}' [% item.dir %]/chr.sizes >> real_chr.csv
[% END -%]

cat chr_length_chrUn.csv real_chr.csv > chr_length.csv
rm real_chr.csv

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
        },
        File::Spec->catfile( $store_dir, "real_chr.sh" )
    ) or die Template->error;

    # id2name.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.name %]
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => [ { name => "S288C", taxon => 4932 }, @data ], },
        File::Spec->catfile( $store_dir, "id2name.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/S288C -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %] -s set01 -p 8 --noaxt -pb lastz --lastz

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "bz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# lpcna
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/lpcna.pl -dt [% data_dir %]/S288C -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "lpcna.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# amp
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/amp.pl -syn -dt [% data_dir %]/S288C -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %] -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => [@data],
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "amp.sh" )
    ) or die Template->error;

    $text = <<'EOF';
cd [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] \
    -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" \
    -a [% data_dir %]/S288Cvs[% item.name %] \
    -lt 10000 -st 0 -ct 0 --parallel 8 --run 1-3,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => [@data],
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "pair_stat.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash

#----------------------------#
# clean RepeatMasker outputs
#----------------------------#
# find [% data_dir %] -name "*.fasta*" | xargs rm

#----------------------------#
# only keeps chr.2bit files
#----------------------------#
# find [% data_dir %] -name "*.fa" | xargs rm

#----------------------------#
# clean pairwise maf
#----------------------------#
find [% data_dir %] -name "mafSynNet" | xargs rm -fr
find [% data_dir %] -name "mafNet" | xargs rm -fr

#----------------------------#
# gzip maf, fas
#----------------------------#
find [% data_dir %] -name "*.maf" | parallel gzip
find [% data_dir %] -name "*.maf.fas" | parallel gzip

#----------------------------#
# clean maf-fasta
#----------------------------#
# rm -fr [% data_dir %]/*_fasta


EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "clean.sh" )
    ) or die Template->error;
}

{    # multi
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast_combine" );
    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {

        # combine
        S288Cvs16_COM => [
            qw{ Spar RM11 YJM789 EC1118 JAY291 Kyokai_no__7 Sigma1278b T7
                DBVPG6044 DBVPG6765 L1528 SK1 UWOPS83 UWOPS87 W303 YPS128 }
        ],
    };

    my @data;
    for my $dbname ( sort keys %{$strains_of} ) {
        my @strains = @{ $strains_of->{$dbname} };
        my $dbs     = join ',', map {"S288Cvs$_"} @strains;
        my $queries = join ',',
            map { $_ . "query" } ( 1 .. scalar @strains - 1 );
        push @data,
            {
            goal_db  => $dbname,
            outgroup => '0query',
            target   => '0target',
            dbs      => $dbs,
            queries  => $queries,
            };
    }

    my $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

[% FOREACH item IN data -%]
# [% item.goal_db %]
perl [% pl_dir %]/alignDB/extra/join_dbs.pl --crude_only \
    --dbs [% item.dbs %] \
    --goal_db [% item.goal_db %] --outgroup [% item.outgroup %] \
    --target [% item.target %] \
    --queries [% item.queries %] \
    --no_insert=1 --trimmed_fasta=1 --length 1000

perl [% pl_dir %]/blastz/refine_fasta.pl \
    --msa mafft -p 8 \
    -i [% data_dir %]/[% item.goal_db %].crude \
    -o [% data_dir %]/[% item.goal_db %]_mafft

perl [% pl_dir %]/tool/catfasta2phyml.pl -f [% data_dir %]/[% item.goal_db %]_mafft/*.fas > [% data_dir %]/all.fasta

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% item.goal_db %] -e yeast_65 -f [% data_dir %]/[% item.goal_db %]_mafft  -lt 1000 -st 0 --parallel 8 --run 1-3,21,40,43

[% END -%]
EOF

    $tt->process(
        \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        File::Spec->catfile( $store_dir, "join.sh" )
    ) or die Template->error;
}

{    # multiz
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast_combine" );
    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {

        S288CvsRM11Spar   => [qw{ Spar RM11 }],
        S288CvsYJM789Spar => [qw{ Spar YJM789 }],
        S288CvsIII        => [qw{ Spar RM11 YJM789 }],

        # from wgs
        S288CvsVIII_WGS =>
            [qw{ Spar RM11 YJM789 EC1118 JAY291 Kyokai_no__7 Sigma1278b T7 }],

        # from sgrp2
        S288CvsXI_SGRP2 => [
            qw{ Spar RM11 YJM789 DBVPG6044 DBVPG6765 L1528 SK1 UWOPS83 UWOPS87
                W303 YPS128 }
        ],

        # combine
        S288CvsXVI_COM => [
            qw{ Spar RM11 YJM789 EC1118 JAY291 Kyokai_no__7 Sigma1278b T7
                DBVPG6044 DBVPG6765 L1528 SK1 UWOPS83 UWOPS87 W303 YPS128 }
        ],

        # branch A
        S288CvsX_A => [
            qw{ Spar RM11 YJM789 EC1118 JAY291 Sigma1278b DBVPG6765 L1528
                UWOPS87 W303 }
        ],

        # branch B
        S288CvsVII_B =>
            [qw{ Spar Kyokai_no__7 T7 DBVPG6044 SK1 UWOPS83 YPS128 }],
    };

    my @data;
    for my $key ( sort keys %{$strains_of} ) {
        my @strains = @{ $strains_of->{$key} };
        push @data,
            {
            out_dir => $key,
            strains => \@strains,
            };
    }

    my $text = <<'EOF';
#!/bin/bash

#----------------------------#
# mz
#----------------------------#
# find [% data_dir %] -name "chrMito*maf" | xargs rm

[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/mz.pl \
    [% FOREACH st IN item.strains -%]
    -d [% data_dir %]/S288Cvs[% st %] \
    [% END -%]
    --tree [% data_dir %]/17way_ml.nwk \
    --out [% data_dir %]/[% item.out_dir %] \
    -syn -p [% parallel %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "mz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#----------------------------#
# maf2fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/maf2fasta.pl \
    --has_outgroup -p [% parallel %] --block \
    -i [% data_dir %]/[% item.out_dir %] \
    -o [% data_dir %]/[% item.out_dir %]_fasta

[% END -%]

#----------------------------#
# mafft
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/refine_fasta.pl \
    --msa mafft --block -p [% parallel %] \
    -i [% data_dir %]/[% item.out_dir %]_fasta \
    -o [% data_dir %]/[% item.out_dir %]_mft

[% END -%]

#----------------------------#
# muscle-quick
#----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/blastz/refine_fasta.pl \
#    --msa muscle --quick --block -p 8 \
#    -i [% data_dir %]/[% item.out_dir %]_fasta \
#    -o [% data_dir %]/[% item.out_dir %]_muscle
#
#[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "maf_fasta.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash

#----------------------------#
# multi_way_batch
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
# mafft
perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl \
    -d [% item.out_dir %] -e yeast_65 \
    --block --id [% data_dir %]/id2name.csv \
    -f [% data_dir %]/[% item.out_dir %]_mft  \
    -lt 5000 -st 0 -ct 0 --parallel [% parallel %] --run common

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "multi.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

if [ ! -d [% data_dir %]/phylo ]
then
    mkdir [% data_dir %]/phylo
fi

#----------------------------#
# concat
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.strains.size > 3 -%]
# [% item.out_dir %]
# concat mafft fas to relaxed phylip
if [ ! -f [% data_dir %]/phylo/[% item.out_dir %].phy ]
then
    perl [% pl_dir %]/blastz/concat_fasta.pl \
        -i [% data_dir %]/[% item.out_dir %]_mft  \
        -o [% data_dir %]/phylo/[% item.out_dir %].phy \
        -p
fi

[% END -%]
[% END -%]

cd [% data_dir %]/phylo

#----------------------------#
# phylo with raxml (ML + rapid bootstrap)
#----------------------------#
[% FOREACH item IN data -%]
[% IF item.strains.size > 3 -%]
# [% item.out_dir %]
if [ -f [% data_dir %]/phylo/[% item.out_dir %].phy.reduced ]
then
    raxml -T 6 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM -O \
        -o Spar -n [% item.out_dir %] \
        -s [% data_dir %]/phylo/[% item.out_dir %].phy.reduced
elif [ -f [% data_dir %]/phylo/[% item.out_dir %].phy ]
then
    raxml -T 6 -f a -m GTRGAMMA -p $RANDOM -N 100 -x $RANDOM \
        -o Spar -n [% item.out_dir %] \
        -s [% data_dir %]/phylo/[% item.out_dir %].phy
fi

[% END -%]
[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "phylo.sh" )
    ) or die Template->error;

}
