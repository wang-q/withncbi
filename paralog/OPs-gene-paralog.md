# Gene with self-aligning information

## Sources

* [gff3 files](https://github.com/wang-q/withncbi/blob/master/ensembl/ensembl.md#gff3)
* [self-aligning](https://github.com/wang-q/withncbi/blob/master/pop/OPs-selfalign.md#arabidopsis)

## Atha

* Convert gff3 to runlists

```bash
GENOME_NAME=Atha

mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}
cd ~/data/alignment/gene-paralog/${GENOME_NAME}

# copy needed files here
cp ~/data/ensembl82/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.29.gff3.gz .
cp ~/data/alignment/self/arabidopsis/Genomes/Atha/chr.sizes chr.sizes
cp ~/data/alignment/self/arabidopsis/Results/Atha/Atha.chr.runlist.yml paralog.yml

# real	50m49.994s
# time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
#     --file Arabidopsis_thaliana.TAIR10.29.gff3.gz \
#     --size chr.sizes \
#     --range 2000 \
#     --clean

# real	0m41.121s
time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
    --file Arabidopsis_thaliana.TAIR10.29.gff3.gz \
    --size chr.sizes \
    --range 2000

# all genes
runlist stat all-gene.yml -s chr.sizes --mk
runlist compare all-gene.yml paralog.yml --op intersect --mk -o all-gene-paralog.yml
runlist stat all-gene-paralog.yml -s chr.sizes --mk

# paralog length
PARALOG_LENGTH=$(runlist stat paralog.yml -s chr.sizes -o stdout | grep 'all,' | cut -d ',' -f 3)

echo "feature,chr,chr_length,feature_size,paralog_length,paralog_size,coverage,paralog_coverage,ratio" \
    > ${GENOME_NAME}.feature.stat.csv
perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t all-gene.yml.csv -m all-gene-paralog.yml.csv \
    -f 0 -f 1 -f2 0 -f2 1 --concat --stdout \
    | cut -d ',' -f 1-4,9 \
    | grep ',all,' \
    | parallel --keep-order echo "{},${PARALOG_LENGTH}" \
    | perl -nl -a -F',' -e '$c1 = $F[3]/$F[2]; $c2 = $F[4]/$F[5]; $r = $c2 / $c1; print join(q{,}, @F[0, 1, 2, 3, 5, 4], $c1, $c2, $r)' \
    >> ${GENOME_NAME}.feature.stat.csv

```

# OsatJap

* Convert gff3 to runlists

```bash
GENOME_NAME=OsatJap

mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}
cd ~/data/alignment/gene-paralog/${GENOME_NAME}

# copy needed files here
cp ~/data/ensembl82/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.29.gff3.gz .
cp ~/data/alignment/self/rice/Genomes/OsatJap/chr.sizes chr.sizes
cp ~/data/alignment/self/rice/Results/OsatJap/OsatJap.chr.runlist.yml paralog.yml

# real	0m40.118s
time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
    --file Oryza_sativa.IRGSP-1.0.29.gff3.gz \
    --size chr.sizes \
    --range 2000

# all genes
runlist stat all-gene.yml -s chr.sizes --mk
runlist compare all-gene.yml paralog.yml --op intersect --mk -o all-gene-paralog.yml
runlist stat all-gene-paralog.yml -s chr.sizes --mk

# paralog length
PARALOG_LENGTH=$(runlist stat paralog.yml -s chr.sizes -o stdout | grep 'all,' | cut -d ',' -f 3)

echo "feature,chr,chr_length,feature_size,paralog_length,paralog_size,coverage,paralog_coverage,ratio" \
    > ${GENOME_NAME}.feature.stat.csv
perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t all-gene.yml.csv -m all-gene-paralog.yml.csv \
    -f 0 -f 1 -f2 0 -f2 1 --concat --stdout \
    | cut -d ',' -f 1-4,9 \
    | grep ',all,' \
    | parallel --keep-order echo "{},${PARALOG_LENGTH}" \
    | perl -nl -a -F',' -e '$c1 = $F[3]/$F[2]; $c2 = $F[4]/$F[5]; $r = $c2 / $c1; print join(q{,}, @F[0, 1, 2, 3, 5, 4], $c1, $c2, $r)' \
    >> ${GENOME_NAME}.feature.stat.csv
```
