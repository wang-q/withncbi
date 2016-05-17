# alt-splice

## Flanking introns

### Put data to `~/data/alt-splice/ase_flanking/`.

```
$ perl ~/Scripts/tool/list_dirtree.pl ~/data/alt-splice/ase_flanking/
+-[ase_flanking]
    |-ASase.pos.txt                               |158522 lines |     3.16M
    |-AScse.pos.txt                               | 49174 lines | 1,003.87K
    |-nonAS.pos.txt                               | 20477 lines |   417.76K
```

### Processing

```bash
cd ~/data/alt-splice/ase_flanking/

# Runtime 20 minutes and 45 seconds.
perl ~/Scripts/alignDB/util/dup_db.pl -f ~/data/dumps/mysql/Human_n12_pop_basic.sql.gz -g Human_n12_ase_flanking

# Runtime 6 hours and 11 minutes.
# ==> Ofg in files: [228173]      Ofg inserted: [212356]
perl ~/Scripts/alignDB/ofg/insert_position.pl \
    -d Human_n12_ase_flanking \
    --style center_intact \
    --parallel 12 \
    --tag intron --type ASase -f ~/data/alt-splice/ase_flanking/ASase.pos.txt \
    --tag intron --type AScse -f ~/data/alt-splice/ase_flanking/AScse.pos.txt \
    --tag intron --type nonAS -f ~/data/alt-splice/ase_flanking/nonAS.pos.txt

perl ~/Scripts/alignDB/init/update_sw_cv.pl \
    -d Human_n12_ase_flanking \
    --parallel 12

perl ~/Scripts/alignDB/init/update_annotation.pl \
    -d Human_n12_ase_flanking \
    -a ~/data/alignment/Ensembl/Human/anno.yml \
    --parallel 12

# Runtime 57 minutes and 44 seconds.
perl ~/Scripts/alignDB/init/update_snp_ld.pl \
    -d Human_n12_ase_flanking --parallel 12

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
    --by type -d Human_n12_ase_flanking \
    --index --chart \
    -o Human_n12_ase_flanking.ofg.xlsx

perl ~/Scripts/alignDB/stat/ld_stat_factory.pl \
    -d Human_n12_ase_flanking \
    --index --chart \
    -o Human_n12_ase_flanking.ld.xlsx

```
