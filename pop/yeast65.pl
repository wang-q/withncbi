#!/usr/bin/perl
use strict;
use warnings;

use Template;

#----------------------------------------------------------#
# RepeatMasker has been done
#----------------------------------------------------------#
my $store_dir = shift
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast65" );

{    # on linux
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast65" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    my $tt   = Template->new;
    my @data = (
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
        {   taxon    => 929587,
            name     => 'CBS_7960',
            coverage => '7.3x 454; 9.69x sanger',
        },
        {   taxon    => 464025,
            name     => 'CLIB215',
            coverage => '6.8x 454; 10.09 sanger',
        },
        {   taxon    => 929629,
            name     => 'CLIB324',
            coverage => '3.2x 454; 3.94x sanger',
        },
        { taxon => 947035, name => 'CLIB382', coverage => '5.96x 454', },
        {   taxon    => 947036,
            name     => 'FL100',
            coverage => '3.2x 454; 3.2x sanger',
        },
        { taxon => 947039, name => 'PW5', coverage => '16.10x 454', },
        { taxon => 929585, name => 'T7',  coverage => '25.4x 454/sanger', },
        { taxon => 471859, name => 'T73', coverage => '13.9x 454', },
        { taxon => 947040, name => 'UC5', coverage => '15.7x 454', },
        {   taxon    => 462210,
            name     => 'Y10',
            coverage => '2.8x 454; 3.81x sanger',
        },
        {   taxon    => 929586,
            name     => 'YJM269',
            coverage => '7.1x 454; 9.59x sanger',
        },

        # wine
        { taxon => 764097, name => 'AWRI796',     coverage => '20x 454', },
        { taxon => 764101, name => 'FostersO',    coverage => '20x 454', },
        { taxon => 764102, name => 'FostersB',    coverage => '20x 454', },
        { taxon => 764098, name => 'Lalvin_QA23', coverage => '20x 454', },
        { taxon => 764099, name => 'Vin13',       coverage => '20x 454', },
        { taxon => 764100, name => 'VL3',         coverage => '20x 454', },

        { taxon => 1095001, name => 'EC9_8', coverage => '30x 454', },
        { taxon => 545124, name => 'AWRI1631', coverage => '7x 454', },
        { taxon => 538975, name => 'M22',      coverage => '2.6x sanger', },
        { taxon => 538976, name => 'YPS163',   coverage => '2.8x sanger', },

        # sgrp data
        { taxon => 900003, name => 'DBVPG6765', coverage => '3x sanger', },
        {   taxon    => 580239,
            name     => 'SK1',
            coverage => '3.27x sanger; 15.61x solexa',
        },
        {   taxon    => 580240,
            name     => 'W303',
            coverage => '2.33x sanger; 3.01x solexa',
        },
        {   taxon    => 900001,
            name     => 'Y55',
            coverage => '3.42x sanger; 8.94x solexa',
        },
    );

    my $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/S288C -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %] -s set01 -p 8 --noaxt -pb lastz --lastz

[% END -%]

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
        File::Spec->catfile( $store_dir, "auto_yeast65_bz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# tar-gzip
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
cd [% data_dir %]/S288Cvs[% item.name %]/

tar -czvf lav.tar.gz   [*.lav   --remove-files
tar -czvf psl.tar.gz   [*.psl   --remove-files
tar -czvf chain.tar.gz [*.chain --remove-files
gzip *.chain
gzip net/*.net
gzip axtNet/*.axt

[% END -%]

#----------------------------#
# clean pairwise maf
#----------------------------#
find [% data_dir %] -name "mafSynNet" | xargs rm -fr
find [% data_dir %] -name "mafNet" | xargs rm -fr

#----------------------------#
# only keeps chr.2bit files
#----------------------------#
# find [% data_dir %] -name "*.fa" | xargs rm
# find [% data_dir %] -name "*.fasta*" | xargs rm

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_yeast65_clean.sh" )
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
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_yeast65_amp.sh" )
    ) or die Template->error;

    $text = <<'EOF';
cd [% data_dir %]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %]
find [% data_dir %]/S288Cvs[% item.name %]/axtNet -name "*.axt.gz" | xargs gzip -d
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/S288Cvs[% item.name %] -lt 1000 -st 1000000 --parallel 8 --run 1-3,21,40
gzip [% data_dir %]/S288Cvs[% item.name %]/axtNet/*.axt

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_yeast65_stat.sh" )
    ) or die Template->error;
}

{    # multi
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast65" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {
        S288CvsALL32 => [
            qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 CBS_7960 CLIB215
                CLIB324 CLIB382 FL100 PW5 T7 T73 UC5 Y10 YJM269 AWRI796 FostersO
                FostersB Lalvin_QA23 Vin13 VL3 EC9_8 Kyokai_no__7 AWRI1631 M22
                YPS163 DBVPG6765 SK1 W303 Y55 }
        ],

        # exclude CLIB382
        S288CvsALL31 => [
            qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 CBS_7960 CLIB215
                CLIB324 FL100 PW5 T7 T73 UC5 Y10 YJM269 AWRI796 FostersO
                FostersB Lalvin_QA23 Vin13 VL3 EC9_8 Kyokai_no__7 AWRI1631 M22
                YPS163 DBVPG6765 SK1 W303 Y55 }
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

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% item.goal_db %] -e yeast_65 -f [% data_dir %]/[% item.goal_db %]_mafft  -lt 1000 -st 1000000 --parallel 8 --run 1-3,21,40

[% END -%]
EOF

    $tt->process(
        \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        File::Spec->catfile( $store_dir, "auto_yeast65_joins.sh" )
    ) or die Template->error;
}

{    # multiz
    my $data_dir = File::Spec->catdir( $ENV{HOME}, "data/alignment/yeast65" );
    my $pl_dir   = File::Spec->catdir( $ENV{HOME}, "Scripts" );

    my $tt         = Template->new;
    my $strains_of = {
        #S288CvsYJM789Spar => [qw{ Spar YJM789 }],
        #S288CvsIII       => [qw{ Spar RM11 YJM789 }],
        #
        ## 10k target length > Spar
        #S288CvsXVIIIGE10m => [
        #    qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 T7 AWRI796 FostersO
        #        FostersB Lalvin_QA23 Vin13 VL3 Kyokai_no__7 DBVPG6765 SK1 W303
        #        Y55 }
        #],

        # 10k target length > Spar,
        # exclude abnormal ins/del ones ( all from wine 454)
        # VL3 Vin13 Lalvin_QA23 AWRI796 FostersO FostersB
        #S288CvsXIIGE10m => [
        #    qw{ Spar AWRI1631 AWRI796 CBS_7960 DBVPG6765 EC1118 EC9_8 FostersB
        #    FostersO JAY291 Kyokai_no__7 Lalvin_QA23 PW5 RM11 SK1 Sigma1278b T7
        #    UC5 VL3 Vin13 W303 Y55 YJM269 YJM789 }
        #],
        #S288CvsXIIGE10m => [
        #    qw{ Spar AWRI1631 CBS_7960 DBVPG6765 EC1118 EC9_8 JAY291
        #    Kyokai_no__7 PW5 RM11 SK1 Sigma1278b T7 UC5 W303 Y55 YJM269 YJM789 }
        #],
        S288CvsVIII => [
            qw{ Spar EC1118 JAY291 Kyokai_no__7 RM11 Sigma1278b T7 YJM789 }
        ],

        ## 1k avg-align > Spar
        #S288CvsXVIIIGE8k => [
        #    qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 PW5 T7 UC5 YJM269
        #        AWRI796 FostersO FostersB Lalvin_QA23 Vin13 VL3 EC9_8
        #        Kyokai_no__7 }
        #],
        #S288CvsXXIIGE8k => [
        #    qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 PW5 T7 UC5 YJM269
        #        AWRI796 FostersO FostersB Lalvin_QA23 Vin13 VL3 EC9_8
        #        Kyokai_no__7 DBVPG6765 SK1 W303 Y55 }
        #],
        #
        ## CLIB382 lacks 4 chr
        #S288CvsXXXI => [
        #    qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 CBS_7960 CLIB215
        #        CLIB324 FL100 PW5 T7 T73 UC5 Y10 YJM269 AWRI796 FostersO
        #        FostersB Lalvin_QA23 Vin13 VL3 EC9_8 Kyokai_no__7 AWRI1631 M22
        #        YPS163 DBVPG6765 SK1 W303 Y55 }
        #],
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
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/mz.pl \
    [% FOREACH st IN item.strains -%]
    -d [% data_dir %]/S288Cvs[% st %] \
    [% END -%]
    --tree [% data_dir %]/33way.nwk \
    --out [% data_dir %]/[% item.out_dir %] \
    -syn -p 8

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_yeast65_mz.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#----------------------------#
# maf2fasta
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/maf2fasta.pl \
    --has_outgroup --id 4932 -p 8 --block \
    -i [% data_dir %]/[% item.out_dir %] \
    -o [% data_dir %]/[% item.out_dir %]_fasta

[% END -%]

#----------------------------#
# mafft
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
perl [% pl_dir %]/blastz/refine_fasta.pl \
    --msa mafft --block -p 8 \
    -i [% data_dir %]/[% item.out_dir %]_fasta \
    -o [% data_dir %]/[% item.out_dir %]_mafft

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

#----------------------------#
# clean
#----------------------------#
[% FOREACH item IN data -%]
# [% item.out_dir %]
cd [% data_dir %]
rm -fr [% item.out_dir %]_fasta

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_yeast65_maf_fasta.sh" )
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
    --block --id 4932 \
    -f [% data_dir %]/[% item.out_dir %]_mafft  \
    -lt 5000 -st 1000000 --parallel 8 --run 1-3,21,40,43

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_yeast65_multi.sh" )
    ) or die Template->error;

}
