#!/usr/bin/perl
use strict;
use warnings;

use Template;

my $data_dir = shift @ARGV;
$data_dir ||= "/home/wangq/data/alignment/yeast";

my $pl_dir = shift @ARGV;
$pl_dir ||= "/home/wangq/Scripts";

{    # wgs
    my $tt   = Template->new;
    my @data = (
        { taxon => 226125, name => 'Spar',   coverage => '7x', },
        { taxon => 285006, name => 'RM11',   coverage => '10x', },
        { taxon => 307796, name => 'YJM789', coverage => '10x', },

        {   taxon    => 574961,
            name     => 'JAY291',
            coverage => '12x 454; 58x solexa s; 95x solexa p',
        },
        {   taxon    => 538975,
            name     => 'Sigma1278b',
            coverage => '45x sanger/solexa',
        },
        { taxon => 643680, name => 'EC1118', coverage => '24x unknown', },

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
        {   taxon    => 721032,
            name     => 'Kyokai_no__7',
            coverage => '9.1x sanger',
        },

        { taxon => 545124, name => 'AWRI1631', coverage => '7x 454', },
        { taxon => 538975, name => 'M22',      coverage => '2.6x sanger', },
        { taxon => 538976, name => 'YPS163',   coverage => '2.8x sanger', },

    );

    my $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# align
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
gunzip -c [% item.name %]/*.gz > [% item.name %]/[% item.name %].fasta
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' [% item.name %]/[% item.name %].fasta
RepeatMasker [% item.name %]/*.fasta -species fungi -xsmall --parallel 4
mv [% item.name %]/[% item.name %].fasta.masked [% item.name %]/[% item.name %].fa
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/S288C_58 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %]_58 -s set11 -p 4
[% END -%]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -e yeast_58 -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/S288Cvs[% item.name %]_58 -lt 5000 -st 1000000 --parallel 4 --run all
[% END -%]

EOF

    $tt->process( \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        "auto_wgs.sh" )
        or die Template->error;
}

{    # sgrp
    my $tt   = Template->new;
    my @data = (
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
# align
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
RepeatMasker [% item.name %]/*.fasta -species fungi -xsmall --parallel 4
mv [% item.name %]/[% item.name %].fasta.masked [% item.name %]/[% item.name %].fa
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/S288C_58 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %]_58 -s set11 -p 4
[% END -%]

#----------------------------#
# stat
#----------------------------#
[% FOREACH item IN data -%]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -e yeast_58 -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/S288Cvs[% item.name %]_58 -lt 5000 -st 1000000 --parallel 4 --run all
[% END -%]

EOF

    $tt->process( \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        "auto_sgrp.sh" )
        or die Template->error;
}

{    # multi
    my $tt         = Template->new;
    my $strains_of = {
        S288CvsYJM789refSpar => [qw{ Spar YJM789 }],
        S288CvsThree         => [qw{ Spar RM11 YJM789 }],
        S288CvsSix           => [qw{ Spar RM11 YJM789 DBVPG6765 SK1 Y55 }],
        S288CvsGE10M18       => [
            qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 T7 AWRI796
                Lalvin_QA23 Vin13 VL3 FostersO FostersB Kyokai_no__7 DBVPG6765
                SK1 Y55 W303
                }
        ],
        S288CvsALL32 => [
            qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 CBS_7960 CLIB215
                CLIB324 FL100 Y10 YJM269 CLIB382 PW5 T7 T73 UC5 AWRI796
                Lalvin_QA23 Vin13 VL3 FostersO FostersB EC9_8 Kyokai_no__7
                AWRI1631 M22 YPS163 DBVPG6765 SK1 Y55 W303
                }
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
perl [% pl_dir %]/alignDB/extra/join_dbs.pl --dbs [% item.dbs %] --goal_db [% item.goal_db %] --outgroup [% item.outgroup %] --target [% item.target %] --queries [% item.queries %] --no_insert=1 --trimmed_fasta=1 --length 1000

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% item.goal_db %] -e yeast_58 -f [% data_dir %]/[% item.goal_db %]  -lt 1000 -st 100000 --parallel 4 --run all

[% END -%]
EOF

    $tt->process( \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        "auto_joins.sh" )
        or die Template->error;
}
