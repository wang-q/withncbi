# alt-splice

## Flanking introns

### Put data to `~/data/alt-splice/ase_flanking/`.

```
$ perl ~/Scripts/tool/list_dirtree.pl ~/data/alt-splice/ase_flanking/
+-[ase_flanking]
    |-ASase.bed                                   |158522 lines |     3.16M
    |-AScse.bed                                   | 49174 lines | 1,003.87K
    |-nonAS.bed                                   | 20477 lines |   417.76K
```

### Processing

```bash
cd ~/data/alt-splice/ase_flanking/

# Runtime 19 minutes and 37 seconds.
perl ~/Scripts/alignDB/util/dup_db.pl -f ~/data/dumps/mysql/Human_n12_pop_basic.sql.gz -g Human_n12_ase_flanking

# Runtime 8 hours and 42 minutes.
# ==> Ofg in files: [228173]      Ofg inserted: [212356]
perl ~/Scripts/alignDB/ofg/insert_bed.pl \
    -d Human_n12_ase_flanking \
    --style center_intact \
    --parallel 12 \
    --tag intron --type ASase -f ~/data/alt-splice/ase_flanking/ASase.bed \
    --tag intron --type AScse -f ~/data/alt-splice/ase_flanking/AScse.bed \
    --tag intron --type nonAS -f ~/data/alt-splice/ase_flanking/nonAS.bed

perl ~/Scripts/alignDB/init/update_sw_cv.pl \
    -d Human_n12_ase_flanking \
    --parallel 12

perl ~/Scripts/alignDB/init/update_feature.pl \
    -d Human_n12_ase_flanking \
    -e homo_sapiens_core_82_37 \
    --parallel 12

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
    --by type -d Human_n12_ase_flanking \
    --index --chart \
    -o Human_n12_ase_flanking.ofg.xlsx

```
