# LD between indels and SNPs

## Yeast

1. dup DB

```bash
perl ~/Scripts/alignDB/util/dup_db.pl -f ~/data/dumps/mysql/Scer_n8_pop.sql.gz -g Scer_n8_pop_LD
```

2. Optional: insert ofg

```bash
perl ~/Scripts/alignDB/ofg/insert_position.pl \
    -d Scer_n8_pop_LD \
    --style center_intact --batch 10 --parallel 8 \
    --tag spo11 --type hotspot -f ~/data/ofg/spo11/spo11_hot.pos.txt

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d Scer_n8_pop_LD --batch 10 --parallel 8
perl ~/Scripts/alignDB/init/update_feature.pl \
    -d Scer_n8_pop_LD \
    -e saccharomyces_cerevisiae_core_29_82_4 \
    --batch 10 --parallel 8
```

3. update ld

```bash
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d Scer_n8_pop_LD --parallel 8 --batch 10
```

4. stat

```bash
perl ~/Scripts/alignDB/stat/ld_stat_factory.pl -d Scer_n8_pop_LD --index --chart
perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl -d Scer_n8_pop_LD --index --chart --run 21
```

## Others

```bash
perl ~/Scripts/alignDB/util/dup_db.pl -d DmelvsXXII -g DmelvsXXII_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d DmelvsXXII_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d AthvsXIX -g AthvsXIX_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d AthvsXIX_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d NipvsXXIV -g NipvsXXIV_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d NipvsXXIV_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d MousevsXIIS -g MousevsXIIS_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d MousevsXIIS_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d HumanvsXI -g HumanvsXI_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d HumanvsXI_LD

perl ld_stat_factory.pl -s 114.212.202.159 -d DmelvsXXII_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d AthvsXIX_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d NipvsXXIV_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d MousevsXIIS_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d HumanvsXI_LD 

# bacs
perl ~/Scripts/alignDB/util/dup_db.pl -d Escherichia_coli -g Escherichia_coli_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d Escherichia_coli_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d Salmonella_enterica -g Salmonella_enterica_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d Salmonella_enterica_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d Staphylococcus_aureus -g Staphylococcus_aureus_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d Staphylococcus_aureus_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d Streptococcus_pneumoniae -g Streptococcus_pneumoniae_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d Streptococcus_pneumoniae_LD

perl ~/Scripts/alignDB/util/dup_db.pl -d Streptococcus_pyogenes -g Streptococcus_pyogenes_LD
perl ~/Scripts/alignDB/init/update_snp_ld.pl -d Streptococcus_pyogenes_LD

perl ld_stat_factory.pl -s 114.212.202.159 -d Escherichia_coli_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d Salmonella_enterica_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d Staphylococcus_aureus_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d Streptococcus_pneumoniae_LD 

perl ld_stat_factory.pl -s 114.212.202.159 -d Streptococcus_pyogenes_LD 

```
