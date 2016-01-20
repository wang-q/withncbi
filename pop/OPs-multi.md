# Multiple alignments

Make sure local ensembl databases have been built.

Create result and dump directories.

```bash
mkdir -p ~/data/alignment/xlsx
mkdir -p ~/data/dumps/mysql
```

## Human (complete genomics data)

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Human_n11cg_chimp_basic \
    -da ~/data/alignment/human_cg/fas_mft  \
    -e homo_sapiens_core_82_37 \
    --block --outgroup \
    -chr ~/data/alignment/human_cg/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run basic

perl ~/Scripts/alignDB/util/dup_db.pl -d Human_n11cg_chimp_basic -f ~/data/dumps/mysql/Human_n11cg_chimp_basic.sql.gz

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Human_n11cg_chimp \
    -da ~/data/alignment/human_cg/fas_mft  \
    -e homo_sapiens_core_82_37 \
    --block --outgroup \
    -chr ~/data/alignment/human_cg/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Human_n11cg_chimp -f ~/data/dumps/mysql/Human_n11cg_chimp.sql.gz

```

## Fly (dpgp82)

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Dmel_n22_pop_basic \
    -da ~/data/alignment/dpgp82/Dmel_n22_pop_refined \
    -e drosophila_melanogaster_core_82_602 \
    --block \
    -chr ~/data/alignment/dpgp82/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run basic

perl ~/Scripts/alignDB/util/dup_db.pl -d Dmel_n22_pop_basic -f ~/data/dumps/mysql/Dmel_n22_pop_basic.sql.gz

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Dmel_n22_pop \
    -da ~/data/alignment/dpgp82/Dmel_n22_pop_refined \
    -e drosophila_melanogaster_core_82_602 \
    --block \
    -chr ~/data/alignment/dpgp82/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Dmel_n22_pop -f ~/data/dumps/mysql/Dmel_n22_pop.sql.gz

```

## Yeast

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Scer_n8_pop \
    -da ~/data/alignment/Fungi/scer_wgs/Scer_n8_pop_refined \
    -e saccharomyces_cerevisiae_core_29_82_4 \
    --block \
    -chr ~/data/alignment/Fungi/scer_wgs/chr_length.csv \
    -lt 5000 --parallel 8 --batch 5 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Scer_n8_pop -f ~/data/dumps/mysql/Scer_n8_pop.sql.gz

```
