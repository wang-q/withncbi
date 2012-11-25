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
    || File::Spec->catdir( $ENV{HOME}, "data/alignment/saccharomyces" );
my $parallel = 12;

{    # on linux
    my $data_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/saccharomyces" );
    my $pl_dir      = File::Spec->catdir( $ENV{HOME}, "Scripts" );
    my $kentbin_dir = File::Spec->catdir( $ENV{HOME}, "bin/x86_64" );

    # ensembl genomes 65
    my $fasta_dir
        = File::Spec->catdir( $ENV{HOME}, "data/alignment/saccharomyces/WGS" );

    my $tt = Template->new;

    my @data = (
        {   taxon    => 1160507,
            name     => "Sarb_H_6",
            sciname  => "Saccharomyces arboricola H-6",
            prefix   => "ALIE01",
            coverage => "50.0x 454; SOLiD",
        },
        {   taxon    => 226231,
            name     => "Sbay_623_6C",
            sciname  => "Saccharomyces bayanus 623-6C",
            prefix   => "AACG02",
            coverage => " ",
        },
        {   taxon    => 226127,
            name     => "Sbay_MCYC_623",
            sciname  => "Saccharomyces bayanus MCYC 623",
            prefix   => "AACA01",
            coverage => " ",
        },
        {   taxon    => 545124,
            name     => "Scer_AWRI1631",
            sciname  => "Saccharomyces cerevisiae AWRI1631",
            prefix   => "ABSV01",
            coverage => " ",
        },
        {   taxon    => 764097,
            name     => "Scer_AWRI796",
            sciname  => "Saccharomyces cerevisiae AWRI796",
            prefix   => "ADVS01",
            coverage => "20x 454",
        },
        {   taxon    => 929587,
            name     => "Scer_CBS_7960",
            sciname  => "Saccharomyces cerevisiae CBS 7960",
            prefix   => "AEWL01",
            coverage => "17.0X 454; ABI 3730",
        },
        {   taxon    => 889517,
            name     => "Scer_CEN_PK113_7D",
            sciname  => "Saccharomyces cerevisiae CEN.PK113-7D",
            prefix   => "AEHG01",
            coverage => "18x 454; llumina",
        },
        {   taxon    => 464025,
            name     => "Scer_CLIB215",
            sciname  => "Saccharomyces cerevisiae CLIB215",
            prefix   => "AEWP01",
            coverage => "16.9X 454; ABI 3730",
        },
        {   taxon    => 929629,
            name     => "Scer_CLIB324",
            sciname  => "Saccharomyces cerevisiae CLIB324",
            prefix   => "AEWM01",
            coverage => "7.14X 454; ABI 3730",
        },
        {   taxon    => 1095001,
            name     => "Scer_EC9_8",
            sciname  => "Saccharomyces cerevisiae EC9-8",
            prefix   => "AGSJ01",
            coverage => "30x 454 Roche GS Junior",
        },
        {   taxon    => 947036,
            name     => "Scer_FL100",
            sciname  => "Saccharomyces cerevisiae FL100",
            prefix   => "AEWO01",
            coverage => "7.1X 454; ABI 3730",
        },
        {   taxon    => 764102,
            name     => "Scer_FostersB",
            sciname  => "Saccharomyces cerevisiae FostersB",
            prefix   => "AEHH01",
            coverage => "20x 454",
        },
        {   taxon    => 764101,
            name     => "Scer_FostersO",
            sciname  => "Saccharomyces cerevisiae FostersO",
            prefix   => "AEEZ01",
            coverage => "20x 454",
        },
        {   taxon   => 574961,
            name    => "Scer_JAY291",
            sciname => "Saccharomyces cerevisiae JAY291",
            prefix  => "ACFL01",
            coverage =>
                "12x 454; 58x Solexa single-end reads; 95x Solexa paired-end 454; Solexa",
        },
        {   taxon    => 721032,
            name     => "Scer_Kyokai_no_7",
            sciname  => "Saccharomyces cerevisiae Kyokai no. 7",
            prefix   => "BABQ01",
            coverage => "9.1x ABI 3730xl",
        },
        {   taxon    => 764098,
            name     => "Scer_Lalvin_QA23",
            sciname  => "Saccharomyces cerevisiae Lalvin QA23",
            prefix   => "ADVV01",
            coverage => "20x 454",
        },
        {   taxon    => 538975,
            name     => "Scer_M22",
            sciname  => "Saccharomyces cerevisiae M22",
            prefix   => "ABPC01",
            coverage => " ",
        },
        {   taxon    => 947039,
            name     => "Scer_PW5",
            sciname  => "Saccharomyces cerevisiae PW5",
            prefix   => "AFDC01",
            coverage => "16.10x 454",
        },
        {   taxon    => 285006,
            name     => "Scer_RM11_1a",
            sciname  => "Saccharomyces cerevisiae RM11-1a",
            prefix   => "AAEG01",
            coverage => " ",
        },
        {   taxon    => 658763,
            name     => "Scer_Sigma1278b",
            sciname  => "Saccharomyces cerevisiae Sigma1278b",
            prefix   => "ACVY01",
            coverage => " ",
        },
        {   taxon    => 929585,
            name     => "Scer_T7",
            sciname  => "Saccharomyces cerevisiae T7",
            prefix   => "AFDE01",
            coverage => "25.4x 454; ABI 3730",
        },
        {   taxon    => 471859,
            name     => "Scer_T73",
            sciname  => "Saccharomyces cerevisiae T73",
            prefix   => "AFDF01",
            coverage => "13.9x 454",
        },
        {   taxon    => 947040,
            name     => "Scer_UC5",
            sciname  => "Saccharomyces cerevisiae UC5",
            prefix   => "AFDD01",
            coverage => "15.7x 454",
        },
        {   taxon    => 764100,
            name     => "Scer_VL3",
            sciname  => "Saccharomyces cerevisiae VL3",
            prefix   => "AEJS01",
            coverage => "20x 454",
        },
        {   taxon    => 764099,
            name     => "Scer_Vin13",
            sciname  => "Saccharomyces cerevisiae Vin13",
            prefix   => "ADXC01",
            coverage => "20x 454",
        },
        {   taxon    => 580240,
            name     => "Scer_W303",
            sciname  => "Saccharomyces cerevisiae W303",
            prefix   => "ALAV01",
            coverage => "40.0x 454",
        },
        {   taxon    => 462210,
            name     => "Scer_Y10",
            sciname  => "Saccharomyces cerevisiae Y10",
            prefix   => "AEWK01",
            coverage => "6.6X 454; ABI 3730",
        },
        {   taxon    => 929586,
            name     => "Scer_YJM269",
            sciname  => "Saccharomyces cerevisiae YJM269",
            prefix   => "AEWN01",
            coverage => "16.7X 454; ABI 3730",
        },
        {   taxon    => 307796,
            name     => "Scer_YJM789",
            sciname  => "Saccharomyces cerevisiae YJM789",
            prefix   => "AAFW02",
            coverage => " ",
        },
        {   taxon    => 1087981,
            name     => "Scer_YJSH1",
            sciname  => "Saccharomyces cerevisiae YJSH1",
            prefix   => "AGAW01",
            coverage => "24x 454 GS FLX Titanium; Sanger",
        },
        {   taxon    => 538976,
            name     => "Scer_YPS163",
            sciname  => "Saccharomyces cerevisiae YPS163",
            prefix   => "ABPD01",
            coverage => " ",
        },
        {   taxon    => 1227742,
            name     => "Scer_ZTW1",
            sciname  => "Saccharomyces cerevisiae ZTW1",
            prefix   => "AMDD01",
            coverage => "20.0x 454; Sanger",
        },
        {   taxon    => 226230,
            name     => "Skud_IFO_1802",
            sciname  => "Saccharomyces kudriavzevii IFO 1802",
            prefix   => "AACI03",
            coverage => "3.4x Sanger",
        },
        {   taxon    => 226126,
            name     => "Smik_IFO_1815",
            sciname  => "Saccharomyces mikatae IFO 1815",
            prefix   => "AACH01",
            coverage => " ",
        },
        {   taxon    => 226125,
            name     => "Spar_NRRL_Y_17217",
            sciname  => "Saccharomyces paradoxus NRRL Y-17217",
            prefix   => "AABY01",
            coverage => " ",
        },
        {   taxon    => 1214527,
            name     => "Spas_CCY48_91",
            sciname  => "Saccharomyces pastorianus CCY48 - 91",
            prefix   => "ALJS01",
            coverage => "30.0x 454 FLX",
        },
        {   taxon    => 520522,
            name     => "Spas_Weihenstephan_34_70",
            sciname  => "Saccharomyces pastorianus Weihenstephan 34/70",
            prefix   => "ABPO01",
            coverage => " ",
        },
    );

    my @files_fasta
        = File::Find::Rule->file->name('*.fasta.gz')->in($fasta_dir);

    for my $item (@data) {

        # match the most similar name
        my ($fasta) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( basename($_), $item->{prefix} ) ] }
            @files_fasta;
        $item->{fasta} = $fasta;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $item->{name} );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    print Dump \@data;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],[% item.sciname.split(' ').slice(0, 1).join(',') %],[% item.name %],,
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
[% item.taxon %],chrUn,999999999,[% item.name %]/WGS
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => [ { name => "S288C", taxon => 4932 }, @data ], },
        File::Spec->catfile( $store_dir, "chr_length_chrUn.csv" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
perl -aln -F"\t" -e 'print qq{[% item.taxon %],$F[0],$F[1],[% item.name %]/WGS}' [% item.dir %]/chr.sizes >> real_chr.csv
[% END -%]

cat chr_length_chrUn.csv real_chr.csv > chr_length.csv
rm real_chr.csv

EOF
    $tt->process(
        \$text,
        {   data     => [ { name => "S288C", taxon => 4932 }, @data ],
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
# unzip, filter and split
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %] 
echo [% item.name %]

if [ ! -d [% item.dir %] ]
then
    mkdir [% item.dir %]
fi

cd [% item.dir %]
gzip -d -c [% item.fasta %] > toplevel.fa
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' toplevel.fa

[% kentbin_dir %]/faCount toplevel.fa > count
cat count | perl -aln -e 'next if $F[0] =~ /^#/; next if $F[0] eq 'total'; print $F[0] if $F[1] > 100000;' | uniq | sort > long_contig
cat count | perl -aln -e 'next if $F[0] =~ /^#/; next if $F[0] eq 'total'; print $F[0] if $F[1] <= 100000 and $F[1] > 5000 and $F[6]/$F[1] < 0.05;' | uniq | sort > short_contig

[% kentbin_dir %]/faSomeRecords toplevel.fa long_contig long_contigs.fa 
[% kentbin_dir %]/faSplit byname long_contigs.fa  .
[% kentbin_dir %]/faSomeRecords toplevel.fa short_contig short_contigs.fa 

rm toplevel.fa count long_contigs.fa long_contig short_contig
rename 's/fa$/fasta/' *.fa;
    
[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "file.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# RepeatMasker
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

cd [% item.dir %]
RepeatMasker [% item.dir %]/*.fasta -species Fungi -xsmall --parallel [% parallel %]

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
        File::Spec->catfile( $store_dir, "rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# RepeatMasker
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
echo [% item.name %]

for i in [% item.dir %]/*.fasta;
do
    if [ -f $i.masked ];
    then
        rename 's/fasta.masked$/fa/' $i.masked;
        find [% item.dir %] -type f -name "`basename $i`*" | xargs rm 
    fi;
done;

[% END -%]

echo Please check the following files
find [% data_dir %] -name "*.fasta"

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir,
            parallel    => $parallel,
        },
        File::Spec->catfile( $store_dir, "clean-rm.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz S288C
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
# S288C versus [% data.$i.name %]
perl [% pl_dir %]/blastz/bz.pl \
    -dt [% data_dir %]/S288C \
    -dq [% data_dir %]/[% data.$i.name %] \
    -dl [% data_dir %]/S288Cvs[% data.$i.name %] \
    -s set01 -p [% parallel %] --noaxt -pb lastz --lastz

[% END -%]

#----------------------------#
# blastz non Scer (Sarb Sbay Skud Smik Spar Spas)
#----------------------------#
[% SET others = [] -%]
[% FOREACH item IN data;
    others.push(item) UNLESS item.name.match('^Scer');
END -%]
[% FOREACH i IN [ 0 .. others.max ] -%]
[% FOREACH j IN [ i .. others.max ] -%]
[% NEXT IF i == j -%]
# [% others.$i.name %] versus [% others.$j.name %]
perl [% pl_dir %]/blastz/bz.pl \
    -dt [% data_dir %]/[% others.$i.name %] \
    -dq [% data_dir %]/[% others.$j.name %] \
    -dl [% data_dir %]/[% others.$i.name %]vs[% others.$j.name %] \
    -s set01 -p [% parallel %] --noaxt -pb lastz --lastz

[% END -%]
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
# lpcna S288C
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
# S288C versus [% data.$i.name %]
perl [% pl_dir %]/blastz/lpcna.pl \
    -dt [% data_dir %]/S288C \
    -dq [% data_dir %]/[% data.$i.name %] \
    -dl [% data_dir %]/S288Cvs[% data.$i.name %] \
    -p [% parallel %] 

[% END -%]

#----------------------------#
# lpcna non Scer (Sarb Sbay Skud Smik Spar Spas)
#----------------------------#
[% SET others = [] -%]
[% FOREACH item IN data;
    others.push(item) UNLESS item.name.match('^Scer');
END -%]
[% FOREACH i IN [ 0 .. others.max ] -%]
[% FOREACH j IN [ i .. others.max ] -%]
[% NEXT IF i == j -%]
# [% others.$i.name %] versus [% others.$j.name %]
perl [% pl_dir %]/blastz/lpcna.pl \
    -dt [% data_dir %]/[% others.$i.name %] \
    -dq [% data_dir %]/[% others.$j.name %] \
    -dl [% data_dir %]/[% others.$i.name %]vs[% others.$j.name %] \
    -p [% parallel %] 

[% END -%]
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
cd [% data_dir %]

#----------------------------#
# amp S288C
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
# S288C versus [% data.$i.name %]
perl [% pl_dir %]/blastz/amp.pl -syn \
    -dt [% data_dir %]/S288C \
    -dq [% data_dir %]/[% data.$i.name %] \
    -dl [% data_dir %]/S288Cvs[% data.$i.name %] \
    -p [% parallel %] 

[% END -%]

#----------------------------#
# amp non Scer (Sarb Sbay Skud Smik Spar Spas)
#----------------------------#
[% SET others = [] -%]
[% FOREACH item IN data;
    others.push(item) UNLESS item.name.match('^Scer');
END -%]
[% FOREACH i IN [ 0 .. others.max ] -%]
[% FOREACH j IN [ i .. others.max ] -%]
[% NEXT IF i == j -%]
# [% others.$i.name %] versus [% others.$j.name %]
perl [% pl_dir %]/blastz/amp.pl -syn \
    -dt [% data_dir %]/[% others.$i.name %] \
    -dq [% data_dir %]/[% others.$j.name %] \
    -dl [% data_dir %]/[% others.$i.name %]vs[% others.$j.name %] \
    -p [% parallel %] 

[% END -%]
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
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# stat S288C
#----------------------------#
[% FOREACH i IN [ 0 .. data.max ] -%]
# S288C versus [% data.$i.name %]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl \
    -d S288Cvs[% data.$i.name %] \
    -t="4932,S288C" \
    -q "[% data.$i.taxon %],[% data.$i.name %]" \
    -a [% data_dir %]/S288Cvs[% data.$i.name %] \
    -at 5000 -st 0 -ct 0 --parallel [% parallel %] --run 1-3,21,40

[% END -%]

#----------------------------#
# stat non Scer (Sarb Sbay Skud Smik Spar Spas)
#----------------------------#
[% SET others = [] -%]
[% FOREACH item IN data;
    others.push(item) UNLESS item.name.match('^Scer');
END -%]
[% FOREACH i IN [ 0 .. others.max ] -%]
[% FOREACH j IN [ i .. others.max ] -%]
[% NEXT IF i == j -%]
# [% others.$i.name %] versus [% others.$j.name %]
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl \
    -d [% others.$i.name %]vs[% others.$j.name %] \
    -t "[% others.$i.taxon %],[% others.$i.name %]" \
    -q "[% others.$j.taxon %],[% others.$j.name %]" \
    -a [% data_dir %]/[% others.$i.name %]vs[% others.$j.name %] \
    -at 5000 -st 0 -ct 0 --parallel [% parallel %] --run 1-3,21,40

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
        File::Spec->catfile( $store_dir, "pair_stat.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash

#----------------------------#
# only keeps chr.2bit files
#----------------------------#
# find [% data_dir %] -name "*.fa" | xargs rm

#----------------------------#
# clean pairwise maf
#----------------------------#
# find [% data_dir %] -name "mafSynNet" | xargs rm -fr
# find [% data_dir %] -name "mafNet" | xargs rm -fr

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
            kentbin_dir => $kentbin_dir,
        },
        File::Spec->catfile( $store_dir, "clean.sh" )
    ) or die Template->error;
}

#{    # multiz
#    my $data_dir
#        = File::Spec->catdir( $ENV{HOME}, "data/alignment/aspergillus" );
#    my $pl_dir = File::Spec->catdir( $ENV{HOME}, "Scripts" );
#
#    my $tt = Template->new;
#    my $strains_of
#        = { AfumvsVII => [qw{ Acla Afla Anid Anig Aory Ater Nfis }], };
#
#    my @data;
#    for my $key ( sort keys %{$strains_of} ) {
#        my @strains = @{ $strains_of->{$key} };
#        push @data,
#            {
#            out_dir => $key,
#            strains => \@strains,
#            };
#    }
#
#    my $text = <<'EOF';
##!/bin/bash
#
##----------------------------#
## mz
##----------------------------#
## find . -name "*MT.synNet*" | xargs rm
#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/blastz/mz.pl \
#    [% FOREACH st IN item.strains -%]
#    -d [% data_dir %]/Afumvs[% st %] \
#    [% END -%]
#    --tree [% data_dir %]/8way.nwk \
#    --out [% data_dir %]/[% item.out_dir %] \
#    -syn -p [% parallel %]
#
#[% END -%]
#
#EOF
#    $tt->process(
#        \$text,
#        {   data     => \@data,
#            data_dir => $data_dir,
#            pl_dir   => $pl_dir,
#            parallel => $parallel,
#        },
#        File::Spec->catfile( $store_dir, "mz.sh" )
#    ) or die Template->error;
#
#    $text = <<'EOF';
##----------------------------#
## maf2fasta
##----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/alignDB/util/maf2fasta.pl \
#    --has_outgroup --id 330879 -p [% parallel %] --block \
#    -i [% data_dir %]/[% item.out_dir %] \
#    -o [% data_dir %]/[% item.out_dir %]_fasta
#
#[% END -%]
#
##----------------------------#
## mafft
##----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
#perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
#    --msa mafft --block -p [% parallel %] \
#    -i [% data_dir %]/[% item.out_dir %]_fasta \
#    -o [% data_dir %]/[% item.out_dir %]_mft
#
#[% END -%]
#
##----------------------------#
## muscle-quick
##----------------------------#
##[% FOREACH item IN data -%]
### [% item.out_dir %]
##perl [% pl_dir %]/alignDB/util/refine_fasta.pl \
##    --msa muscle --quick --block -p [% parallel %] \
##    -i [% data_dir %]/[% item.out_dir %]_fasta \
##    -o [% data_dir %]/[% item.out_dir %]_mslq
##
##[% END -%]
#
#EOF
#
#    $tt->process(
#        \$text,
#        {   data     => \@data,
#            data_dir => $data_dir,
#            pl_dir   => $pl_dir,
#            parallel => $parallel,
#        },
#        File::Spec->catfile( $store_dir, "maf_fasta.sh" )
#    ) or die Template->error;
#
#    $text = <<'EOF';
##!/bin/bash
#
##----------------------------#
## multi_way_batch
##----------------------------#
#[% FOREACH item IN data -%]
## [% item.out_dir %]
## mafft
#perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl \
#    -d [% item.out_dir %] -e Afum_65 \
#    --block --id 330879 \
#    -f [% data_dir %]/[% item.out_dir %]_mft  \
#    -lt 5000 -st 0 -ct 0 --parallel [% parallel %] --run 1-3,21,40
#
#[% END -%]
#
#EOF
#    $tt->process(
#        \$text,
#        {   data     => \@data,
#            data_dir => $data_dir,
#            pl_dir   => $pl_dir,
#            parallel => $parallel,
#        },
#        File::Spec->catfile( $store_dir, "multi.sh" )
#    ) or die Template->error;
#}
