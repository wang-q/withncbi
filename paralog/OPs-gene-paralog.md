# Gene with self-aligning information

## Sources

* [gff3 files](https://github.com/wang-q/withncbi/blob/master/ensembl/ensembl.md#gff3)
* [self-aligning](https://github.com/wang-q/withncbi/blob/master/pop/OPs-selfalign.md#arabidopsis)

## MITE

* http://pmite.hzau.edu.cn/download_mite/
* Downloaded files in `~/data/alignment/gene-paralog/download_mite/MITE/MITE-SEQ-V2/`.

```bash
mkdir -p ~/data/alignment/gene-paralog/
cd ~/data/alignment/gene-paralog/
perl ~/Scripts/download/list.pl -ncp -u http://pmite.hzau.edu.cn/download_mite/
perl ~/Scripts/download/download.pl -i download_mite.yml
```

## Scripts

### `stat_runlists.sh`

Need three arguments:
* Multikey runlist
* A genomic feature, such as `paralog.yml`
* chr.sizes

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/stat_runlists.sh
#!/bin/bash

FILE_Y1="$1"
FILE_Y2="$2"
FILE_SIZE="$3"

BASE_Y1=`basename "${FILE_Y1%.*}"`
BASE_Y2=`basename "${FILE_Y2%.*}"`
FILE_OUTPUT="stat.${BASE_Y1}.${BASE_Y2}.csv"

echo "==> parameters"
echo "    " $@

echo "==> compare ${BASE_Y1} ${BASE_Y2}"
runlist stat2 ${FILE_Y1} ${FILE_Y2} -s ${FILE_SIZE} --op intersect --mk -o stdout \
    > ${FILE_OUTPUT}.tmp

echo "==> output ${FILE_OUTPUT}"
head -n 1 ${FILE_OUTPUT}.tmp \
    | cut -d ',' -f 1,3-9 \
    > ${FILE_OUTPUT}
cat ${FILE_OUTPUT}.tmp \
    | grep ',all,' \
    | cut -d ',' -f 1,3-9 \
    >> ${FILE_OUTPUT}

echo "==> clean"
rm ${FILE_OUTPUT}.tmp

EOF

chmod +x ~/data/alignment/gene-paralog/stat_runlists.sh

```

### `concat_csv.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/concat_csv.sh
#!/bin/bash

USAGE="Need two or more csv files.
Usage: $0 csv1 csv2 ... csvN"

if [ "$#" -lt 2 ]; then
	echo "$USAGE"
	exit 1
fi

echo "==> parameters"
echo "    " $@

# First file
COUNTER=1
cp "$1" tmp.concat.${COUNTER}.csv
shift

# Additional files
while (( "$#" ));
do
    COUNTER_1=${COUNTER}
    let "COUNTER++"  
    COUNTER_2=${COUNTER}
    perl ~/Scripts/alignDB/util/merge_csv.pl \
        -t tmp.concat.${COUNTER_1}.csv -m $1 \
        -f 0 -f2 0 --concat --stdout \
        > tmp.concat.${COUNTER_2}.csv
    shift
done

cp tmp.concat.${COUNTER}.csv concat.csv
rm tmp.concat.*.csv

EOF

chmod +x ~/data/alignment/gene-paralog/concat_csv.sh

```

### `proc_all_gene.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/proc_all_gene.sh
#!/bin/bash

USAGE="Usage: $0 GENOME_NAME FEATURE"

if [ "$#" -lt 2 ]; then
	echo "$USAGE"
	exit 1
fi

echo "====> parameters"
echo "    " $@

GENOME_NAME=$1
FEATURE=$2

cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

../../stat_runlists.sh ../feature/all-gene.yml ../data/${FEATURE}.yml ../data/chr.sizes
mv stat.all-gene.${FEATURE}.csv ${GENOME_NAME}.all-gene.${FEATURE}.csv

EOF

```

### `proc_sep_gene.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/proc_sep_gene.sh
#!/bin/bash

USAGE="Usage: $0 GENOME_NAME FEATURE"

if [ "$#" -lt 2 ]; then
	echo "$USAGE"
	exit 1
fi

echo "====> parameters"
echo "    " $@

GENOME_NAME=$1
FEATURE=$2

cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

for ftr in gene upstream downstream exon five_prime_UTR three_prime_UTR CDS intron
do
    echo "====> ${ftr}"
    sleep 1
    ../../stat_runlists.sh ../feature/sep-${ftr}.yml ../data/${FEATURE}.yml ../data/chr.sizes
    cat stat.sep-${ftr}.${FEATURE}.csv \
        | cut -d ',' -f 1,3,5 \
        > stat.sep-${ftr}.${FEATURE}.csv.tmp
done

echo "====> concat gene"
../../concat_csv.sh stat.sep-gene.${FEATURE}.csv.tmp stat.sep-upstream.${FEATURE}.csv.tmp stat.sep-downstream.${FEATURE}.csv.tmp
echo "gene_id,gene_length,gene_${FEATURE},upstream,upstream_${FEATURE},downstream,downstream_${FEATURE}" \
    > ${GENOME_NAME}.gene.${FEATURE}.csv
cat concat.csv \
    | cut -d ',' -f 1,2,3,5,6,8,9 \
    >> ${GENOME_NAME}.gene.${FEATURE}.csv
rm concat.csv

echo "====> concat trans"
../../concat_csv.sh stat.sep-exon.${FEATURE}.csv.tmp stat.sep-intron.${FEATURE}.csv.tmp stat.sep-CDS.${FEATURE}.csv.tmp stat.sep-five_prime_UTR.${FEATURE}.csv.tmp stat.sep-three_prime_UTR.${FEATURE}.csv.tmp
echo "trans_id,exon_length,exon_${FEATURE},intron,intron_${FEATURE},CDS,CDS_${FEATURE},five_prime_UTR,five_prime_UTR_${FEATURE},three_prime_UTR,three_prime_UTR_${FEATURE}" \
    > ${GENOME_NAME}.trans.${FEATURE}.csv
cat concat.csv \
    | cut -d ',' -f 1,2,3,5,6,8,9,11,12,14,15,17,18 \
    >> ${GENOME_NAME}.trans.${FEATURE}.csv
rm concat.csv

echo "====> clean"
rm stat.sep-*.${FEATURE}.csv.tmp
rm stat.sep-*.${FEATURE}.csv

EOF

```

### `proc_mite.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/proc_mite.sh
#!/bin/bash

USAGE="Usage: $0 GENOME_NAME"

if [ "$#" -lt 1 ]; then
	echo "$USAGE"
	exit 1
fi

echo "====> parameters"
echo "    " $@

GENOME_NAME=$1
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

echo "==> mite_stat"
faops size mite.fa \
    | perl -nla -F"\t" -MStatistics::Descriptive -e '
        BEGIN {$stat = Statistics::Descriptive::Full->new;}
        next unless defined $F[1];
        $stat->add_data($F[1]);
        END {
            printf qq{Total:\t%d\nSum:\t%d\n}, $stat->count, $stat->sum;
            printf qq{Median:\t%d\nMean:\t%.1f\n}, $stat->median, $stat->mean;
            printf qq{Min:\t%d\nMax:\t%d\n}, $stat->min, $stat->max;
            printf qq{\n};
            printf qq{Length\t#Seqs\n};
            %distrib = $stat->frequency_distribution(10);
            for (sort {$a <=> $b} keys %distrib) {
            	printf qq{%.1f\t%d\n}, $_, $distrib{$_};
            }
            printf qq{\n};
        }' \
    > mite_stat.txt

echo "==> genome blast"
faops filter -a 40  mite.fa mite.filter.fa
perl ~/Scripts/egas/fasta_blastn.pl -f mite.filter.fa -g genome.fa -o mite.bg.blast --parallel 8
perl ~/Scripts/egas/blastn_genome.pl -f mite.bg.blast -g genome.fa -o mite.bg.fasta -c 0.95 --parallel 8
cat mite.fa mite.bg.fasta \
    | faops filter -u stdin stdout \
    | faops filter -a 40 stdin stdout \
    > mite.all.fasta

echo "==> sparsemem_exact"
perl ~/Scripts/egas/sparsemem_exact.pl -f mite.all.fasta -l 40 -g genome.fa -o mite.replace.tsv
cat mite.replace.tsv \
    | perl -nla -F"\t" -e ' print for @F' \
    | grep ':' \
    | sort | uniq \
    > mite.position.txt

echo "==> mite covers"
wc -l mite.position.txt >> mite_stat.txt
runlist covers mite.position.txt -o mite.yml
runlist stat mite.yml -s chr.sizes

echo "==> clean"
rm mite.all.fasta mite.bg.blast mite.bg.fasta mite.filter.fa mite.replace.tsv

EOF

```

## Atha

1. Prepare

    ```bash
    GENOME_NAME=Atha

    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

    # copy needed files here
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    cp ~/data/ensembl82/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.29.gff3.gz gff3.gz
    cp ~/data/alignment/self/arabidopsis/Genomes/Atha/chr.sizes chr.sizes
    cp ~/data/alignment/self/arabidopsis/Results/Atha/Atha.chr.runlist.yml paralog.yml
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/03_arabidopsis_mite_seq.fa -O mite.fa

    # genome
    find ~/data/alignment/Ensembl/${GENOME_NAME} -type f -name "*.fa" \
        | sort | xargs cat \
        | perl -nl -e '/^>/ or $_ = uc; print' \
        > genome.fa

    runlist stat paralog.yml -s chr.sizes

    # Convert gff3 to runlists
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/feature

    # real	0m41.121s
    time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
        --file ../data/gff3.gz \
        --size ../data/chr.sizes \
        --range 2000

    # real	50m49.994s
    # time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
    #     --file Arabidopsis_thaliana.TAIR10.29.gff3.gz \
    #     --size chr.sizes \
    #     --range 2000 \
    #     --clean
    ```

2. Gene-paralog stats

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha paralog

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha paralog
    ```

3. MITE

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/data

    bash ~/data/alignment/gene-paralog/proc_mite.sh Atha
    ```

    Don't run in parallel with step 2.

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha mite

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha mite
    ```

## OsatJap

1. Prepare

    ```bash
    GENOME_NAME=OsatJap

    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

    # copy needed files here
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    cp ~/data/ensembl82/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.29.gff3.gz gff3.gz
    cp ~/data/alignment/self/rice/Genomes/OsatJap/chr.sizes chr.sizes
    cp ~/data/alignment/self/rice/Results/OsatJap/OsatJap.chr.runlist.yml paralog.yml
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/26_nipponbare_mite_seq.fa -O mite.fa

    # genome
    find ~/data/alignment/Ensembl/${GENOME_NAME} -type f -name "*.fa" \
        | sort | xargs cat \
        | perl -nl -e '/^>/ or $_ = uc; print' \
        > genome.fa

    # Convert gff3 to runlists
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/feature

    # real	1m15.666s
    time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
        --file ../data/gff3.gz \
        --size ../data/chr.sizes \
        --range 2000
    ```

2. Gene-paralog stats

    ```bash
    cd ~/data/alignment/gene-paralog/OsatJap/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh OsatJap paralog

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh OsatJap paralog
    ```

3. MITE

    ```bash
    cd ~/data/alignment/gene-paralog/OsatJap/data

    bash ~/data/alignment/gene-paralog/proc_mite.sh OsatJap
    ```

    ```bash
    cd ~/data/alignment/gene-paralog/OsatJap/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh OsatJap mite

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh OsatJap mite
    ```
