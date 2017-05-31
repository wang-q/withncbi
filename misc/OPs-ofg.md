# Other features of genome (OFG)

[TOC levels=1-3]: # " "
- [Other features of genome (OFG)](#other-features-of-genome-ofg)
- [spo11 (cell 2011)](#spo11-cell-2011)
    - [spo11: data source](#spo11-data-source)
    - [S288Cvsself center_intact](#s288cvsself-center-intact)
    - [Scer_n8_pop center_intact](#scer-n8-pop-center-intact)
    - [S288Cvsself edge](#s288cvsself-edge)
- [S288c X YJM789 (nature 2008)](#s288c-x-yjm789-nature-2008)
    - [Get data](#get-data)
    - [S288Cvsself center](#s288cvsself-center)
    - [nat08 style edge S288CvsYJM789](#nat08-style-edge-s288cvsyjm789)
    - [nat08 style center S288CvsYJM789](#nat08-style-center-s288cvsyjm789)


# spo11 (cell 2011)

## spo11: data source

Original paper:

Pan, J. et al. A Hierarchical Combination of Factors Shapes the Genome-wide Topography of Yeast
Meiotic Recombination Initiation. Cell 144, 719-731 (2011).

```bash
mkdir -p ~/data/ofg/spo11
cd ~/data/ofg/spo11

wget -N http://www.cell.com/cms/attachment/612087/4902311/mmc2.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f mmc2.xls --sheet 'Spo11 Hot Spot Annotations' \
    | perl -F/,/ -an -e '
    /^CHROM/ and next;
    /^chr/ or next;
    $F[0] =~ s/^chr//;
    if ($F[1] == $F[2]) {
        print qq{$F[0]:$F[1]\n};
    }
    else {
        print qq{$F[0]:$F[1]-$F[2]\n};
    }
    ' \
    > spo11_hot.pos.txt
```

## S288Cvsself center_intact

```bash
cd ~/data/ofg/spo11

perl ~/Scripts/alignDB/util/dup_db.pl -g S288cvsself_spo11 -f ~/data/dumps/mysql/S288cvsself.sql.gz

perl ~/Scripts/alignDB/init/insert_position.pl \
    -d S288cvsself_spo11 --batch 50 --parallel 8 \
    --style center_intact \
    --tag spo11 --type hotspot -f ~/data/ofg/spo11/spo11_hot.pos.txt

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d S288cvsself_spo11 --batch 1 --parallel 8

perl ~/Scripts/alignDB/init/update_annotation.pl \
    -d S288cvsself_spo11 \
    -a ~/data/alignment/Ensembl/S288c/anno.yml \
    --batch 10 --parallel 8

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
    --by tt --index --chart \
    -d S288cvsself_spo11

```

## Scer_n8_pop center_intact

```bash
cd ~/data/ofg/spo11

perl ~/Scripts/alignDB/util/dup_db.pl -g Scer_n7_pop_spo11 -f ~/data/dumps/mysql/Scer_n7_pop.sql.gz

perl ~/Scripts/alignDB/init/insert_position.pl \
    -d Scer_n7_pop_spo11 --batch 10 --parallel 8 \
    --style center_intact \
    --tag spo11 --type hotspot -f ~/data/ofg/spo11/spo11_hot.pos.txt

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d Scer_n7_pop_spo11 --batch 10 --parallel 8

perl ~/Scripts/alignDB/init/update_annotation.pl \
    -d Scer_n7_pop_spo11 \
    -a ~/data/alignment/Ensembl/S288c/anno.yml \
    --batch 10 --parallel 8

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
    --by tt --index --chart \
    -d Scer_n7_pop_spo11
```

## S288Cvsself edge

```bash
cd ~/data/ofg/spo11

perl ~/Scripts/alignDB/util/dup_db.pl -g S288cvsself_spo11_edge -f ~/data/dumps/mysql/S288cvsself.sql.gz

perl ~/Scripts/alignDB/init/insert_position.pl \
    -d S288cvsself_spo11_edge --batch 50 --parallel 8 \
    --style edge \
    --tag spo11 --type hotspot -f ~/data/ofg/spo11/spo11_hot.pos.txt

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d S288cvsself_spo11_edge --batch 1 --parallel 8

perl ~/Scripts/alignDB/init/update_annotation.pl \
    -d S288cvsself_spo11_edge \
    -a ~/data/alignment/Ensembl/S288c/anno.yml \
    --batch 10 --parallel 8

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
    --by tt --index --chart \
    -d S288cvsself_spo11_edge
```

# S288c X YJM789 (nature 2008)

## Get data

Original paper:

Mancera, E., Bourgon, R., Brozzi, A., Huber, W. & Steinmetz, L. M. High-resolution mapping of
meiotic crossovers and non-crossovers in yeast. Nature 454, 479-485 (2008).

```bash
mkdir -p ~/data/ofg/nat08
cd ~/data/ofg/nat08

wget -N http://www.nature.com/nature/journal/v454/n7203/extref/nature07135-s2.zip

unzip nature07135-s2.zip

rm gce_*.bed

perl -anl -F"\t" -e '
    BEGIN { use Path::Tiny; use Roman; }
    $F[0] =~ /chr(\d+)/ or next;
    $chr = Roman($1);
    $F[1] =~ s/\.\d+//;
    $F[2] =~ s/\.\d+//;
    $str = join(qq{\t}, $chr, $F[1], $F[2]);
    path(q{gce_} . $F[6] . q{.bed})->append($str);
    ' nature07135-s2/event_intervals.txt
```

## S288Cvsself center

```bash
cd ~/data/ofg/nat08

perl ~/Scripts/alignDB/util/dup_db.pl -g S288cvsself_nat08 -f ~/data/dumps/mysql/S288cvsself.sql.gz

perl ~/Scripts/alignDB/ofg/insert_bed.pl \
    -d S288cvsself_nat08 \
    --style center --batch 1 --parallel 8 \
    --dG \
    --tag gce --type C  -f ~/data/ofg/nat08/gce_C.bed  \
    --tag gce --type CC -f ~/data/ofg/nat08/gce_CC.bed \
    --tag gce --type M  -f ~/data/ofg/nat08/gce_M.bed  \
    --tag gce --type X  -f ~/data/ofg/nat08/gce_X.bed

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d S288cvsself_nat08 --batch 1 --parallel 8
perl ~/Scripts/alignDB/init/update_feature.pl \
    -d S288cvsself_nat08 \
    -e saccharomyces_cerevisiae_core_29_82_4 \
    --batch 1 --parallel 8

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl -d S288cvsself_nat08 --by tt --index --chart
```

## nat08 style edge S288CvsYJM789

## nat08 style center S288CvsYJM789
