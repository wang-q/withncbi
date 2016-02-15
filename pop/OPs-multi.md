# Multiple alignments

Make sure local ensembl databases have been built.

Create result and dump directories.

```bash
mkdir -p ~/data/alignment/xlsx
mkdir -p ~/data/dumps/mysql
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

## Arabidopsis

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Ath_n19_pop \
    -da ~/data/alignment/arabidopsis82/Ath_n19_pop_refined \
    -e arabidopsis_thaliana_core_29_82_10 \
    --block \
    -chr ~/data/alignment/arabidopsis82/chr_length.csv \
    -lt 5000 --parallel 8 --batch 5 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Ath_n19_pop -f ~/data/dumps/mysql/Ath_n19_pop.sql.gz

```

## Rice

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d OsatJap_n24_pop \
    -da ~/data/alignment/rice82/OsatJap_n24_pop_refined \
    -e oryza_sativa_core_29_82_7 \
    --block \
    -chr ~/data/alignment/rice82/chr_length.csv \
    -lt 5000 --parallel 8 --batch 5 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d OsatJap_n24_pop -f ~/data/dumps/mysql/OsatJap_n24_pop.sql.gz

```

## Fly (dpgp82)

```bash
cd ~/data/alignment/xlsx

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

## Cele

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Cele_n41_pop \
    -da ~/data/alignment/cele82/Cele_n41_pop_refined \
    -e caenorhabditis_elegans_core_29_82_245 \
    --block \
    -chr ~/data/alignment/cele82/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Cele_n41_pop -f ~/data/dumps/mysql/Cele_n41_pop.sql.gz
```

## Ddis

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Ddis_n18_pop \
    -da ~/data/alignment/Protists/Ddis/Ddis_n18_pop_refined \
    -e dictyostelium_discoideum_core_29_82_1 \
    --block \
    --chr ~/data/alignment/Protists/Ddis/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Ddis_n18_pop -f ~/data/dumps/mysql/Ddis_n18_pop.sql.gz
```

## Human (complete genomics data)

```bash
cd ~/data/alignment/xlsx

# Runtime 12 hours and 5 minutes.
perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Human_n11cg_chimp_basic \
    -da ~/data/alignment/human_cg/fas_mft  \
    -e homo_sapiens_core_82_37 \
    --block --outgroup \
    -chr ~/data/alignment/human_cg/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run basic

# Runtime 58 minutes and 43 seconds.
perl ~/Scripts/alignDB/util/dup_db.pl -d Human_n11cg_chimp_basic -f ~/data/dumps/mysql/Human_n11cg_chimp_basic.sql.gz

# Runtime 2 days and 14 hours.
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

## Mouse

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Mouse_n11_pop \
    -da ~/data/alignment/mouse82/Mouse_n11_pop_refined  \
    -e mus_musculus_core_82_38 \
    --block \
    --chr ~/data/alignment/mouse82/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Mouse_n11_pop -f ~/data/dumps/mysql/Mouse_n11_pop.sql.gz

```

## Pfal

```bash
cd ~/data/alignment/xlsx

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d Pfal_n7_pop \
    -da ~/data/alignment/Protists/pfal/Pfal_n7_pop_refined  \
    -e plasmodium_falciparum_core_29_82_3 \
    --block \
    --chr ~/data/alignment/Protists/pfal/chr_length.csv \
    -lt 5000 --parallel 12 \
    --run all

perl ~/Scripts/alignDB/util/dup_db.pl -d Pfal_n7_pop -f ~/data/dumps/mysql/Pfal_n7_pop.sql.gz

```
