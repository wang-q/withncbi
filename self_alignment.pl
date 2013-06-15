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

my $store_dir = shift
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/self_alignment" );
my $parallel = 12;

{    # on linux
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/self_alignment" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    my @data = (
        {   taxon   => 3702,
            name    => "ath",
            ensembl => "ath_65",
            join    => "AthvsLyrata",
        },
        {   taxon   => 7227,
            name    => "Dmel",
            ensembl => "fly_65",
            join    => "DmelvsDsim",
        },
        {   taxon   => 9606,
            name    => "human",
            ensembl => "human_65",
            join    => "HumanvsChimp",
        },
        {   taxon   => 10090,
            name    => "mouse",
            ensembl => "mouse_65",
            join    => "MousevsSpretus_Ei",
        },
        {   taxon   => 39947,
            name    => "nip",
            ensembl => "nip_65",
            join    => "Nipvs9311",
        },
        {   taxon   => 4932,
            name    => "S288C",
            ensembl => "yeast_65",
            join    => "S288CvsSpar",
        },
    );

    for my $item ( sort @data ) {
        my $dir = File::Spec->catdir( $data_dir, $item->{name} );
        $item->{dir} = $dir;
    }

    my $tt = Template->new;
    my $text;

    # chr_length_chrUn.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]/ensembl65
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
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]/ensembl65}' [% item.dir %]/chr.sizes >> real_chr.csv
[% END -%]

cat chr_length_chrUn.csv real_chr.csv > chr_length.csv
rm real_chr.csv

echo Run the following cmds to merge csv files
echo
echo perl [% pl_dir %]/alignDB/util/merge_csv.pl -t [% pl_dir %]/alignDB/init/taxon.csv -m [% data_dir %]/taxon.csv
echo
echo perl [% pl_dir %]/alignDB/util/merge_csv.pl -t [% pl_dir %]/alignDB/init/chr_length.csv -m [% data_dir %]/chr_length.csv
echo

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
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
        { data => \@data, },
        File::Spec->catfile( $store_dir, "id2name.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %]
perl [% pl_dir %]/blastz/bz.pl --parallel [% parallel %] \
    --is_self \
    -s set01 -C 0 --noaxt -pb lastz --lastz \
    -dt [% item.dir %] \
    -dq [% item.dir %] \
    -dl [% data_dir %]/[% item.name %]vsselfalign

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
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
# [% item.name %] 
perl [% pl_dir %]/blastz/lpcna.pl --parallel [% parallel %] \
    -dt [% item.dir %] \
    -dq [% item.dir %] \
    -dl [% data_dir %]/[% item.name %]vsselfalign

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "lpcna.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
    
#----------------------------#
# amp
#----------------------------#
# Don't use -syn here

[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
perl [% pl_dir %]/blastz/amp.pl --parallel [% parallel %] \
    -dt [% item.dir %] \
    -dq [% item.dir %] \
    -dl [% data_dir %]/[% item.name %]vsselfalign

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
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
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl --parallel [% parallel %] \
    -t "[% item.taxon %],[% item.name %]" \
    -q "[% item.taxon %],[% item.name %]" \
    -d [% item.name %]vsselfalign \
    --dir [% data_dir %]/[% item.name %]vsselfalign \
    -e [% item.ensembl %] \
    -lt 1000 --run gc

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "pair_stat.sh" )
    ) or die Template->error;
    
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

[% FOREACH item IN data -%]
# [% item.goal_db %]
perl [% pl_dir %]/alignDB/extra/join_dbs.pl --crude_only \
    --dbs [% item.name %]vsselfalign,[% item.join %] \
    --goal_db [% item.name %]vsselfalign_[% item.join %] \
    --target 0target \
    --queries 0query \
    --outgroup 1query \
    --block --no_insert --trimmed_fasta --length 1000

perl [% pl_dir %]/blastz/refine_fasta.pl --parallel [% parallel %] \
    --outgroup --block \
    --msa mafft \
    -i [% data_dir %]/[% item.name %]vsselfalign_[% item.join %].crude \
    -o [% data_dir %]/[% item.name %]vsselfalign_[% item.join %]_mft

#perl [% pl_dir %]/tool/catfasta2phyml.pl -f [% data_dir %]/[% item.goal_db %]_mafft/*.fas > [% data_dir %]/all.fasta

#perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% item.goal_db %] -e yeast_65 -f [% data_dir %]/[% item.goal_db %]_mafft  -lt 1000 -st 0 --parallel 8 --run 1-3,21,40,43

[% END -%]
EOF

    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
            parallel => $parallel,
        },
        File::Spec->catfile( $store_dir, "join.sh" )
    ) or die Template->error;
}

