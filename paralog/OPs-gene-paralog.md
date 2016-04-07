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

echo "==> bases"
echo "    " ${BASE_Y1} ${BASE_Y2}

echo "==> stat ${BASE_Y2}"
LENGTH_Y2=$(runlist stat ${FILE_Y2} -s ${FILE_SIZE} -o stdout | grep 'all,' | cut -d ',' -f 3)

echo "==> stat ${BASE_Y1}"
runlist stat ${FILE_Y1} -s ${FILE_SIZE} --mk -o stdout \
    | grep ',all,' \
    | cut -d ',' -f 1,3,4 \
    | parallel -j 1 --keep-order echo "{},${LENGTH_Y2}" \
    > tmp1.csv

echo "==> compare ${BASE_Y2}"
runlist compare ${FILE_Y1} ${FILE_Y2} --op intersect --mk -o tmp_intersect.yml

echo "==> stat ${BASE_Y1} ${BASE_Y2}"
runlist stat tmp_intersect.yml -s ${FILE_SIZE} --mk -o stdout \
    | grep ',all,' \
    | cut -d ',' -f 1,4 \
    > tmp2.csv

echo "==> concat csv"
echo "key,chr_length,key_size,${BASE_Y2}_length,key,${BASE_Y2}_size," \
    > tmp_output.csv
perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t tmp1.csv -m tmp2.csv \
    -f 0 -f2 0 --concat --stdout \
    >> tmp_output.csv

echo "==> output"
echo "    " ${FILE_OUTPUT}
cat tmp_output.csv \
    | cut -d ',' -f 1-4,6 \
    > ${FILE_OUTPUT}

echo "==> clean"
rm tmp_intersect.yml tmp1.csv tmp2.csv tmp_output.csv

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
mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
FEATURE=$2

../../stat_runlists.sh ../feature/all-gene.yml ../data/${FEATURE}.yml ../data/chr.sizes
cat stat.all-gene.${FEATURE}.csv \
    | perl -nla -F',' -e '
        ($in_rows) = grep {/^\d+$/} @F;
        if ($in_rows) {
            $c1 = $F[2]/$F[1];
            $c2 = $F[4]/$F[3];
            $r = $c2 / $c1;
            print join(q{,}, @F, $c1, $c2, $r)
        }
        else {
            print join(q{,}, @F, qw{c1 c2 ratio})
        }
        ' \
    > ${GENOME_NAME}.feature.${FEATURE}.csv

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
mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
FEATURE=$2

for ftr in gene upstream downstream exon five_prime_UTR three_prime_UTR CDS intron
do
    echo "====> ${ftr}"
    ../../stat_runlists.sh ../feature/sep-${ftr}.yml ../data/${FEATURE}.yml ../data/chr.sizes
done

../../concat_csv.sh stat.sep-gene.${FEATURE}.csv stat.sep-upstream.${FEATURE}.csv stat.sep-downstream.${FEATURE}.csv
echo "gene_id,gene_length,gene_${FEATURE},upstream,upstream_${FEATURE},downstream,downstream_${FEATURE}" \
    > ${GENOME_NAME}.gene.${FEATURE}.csv
cat concat.csv \
    | cut -d ',' -f 1,3,5,8,10,13,15 \
    >> ${GENOME_NAME}.gene.${FEATURE}.csv
rm concat.csv

../../concat_csv.sh stat.sep-exon.${FEATURE}.csv stat.sep-intron.${FEATURE}.csv stat.sep-CDS.${FEATURE}.csv stat.sep-five_prime_UTR.${FEATURE}.csv stat.sep-three_prime_UTR.${FEATURE}.csv
echo "trans_id,exon_length,exon_${FEATURE},intron,intron_${FEATURE},CDS,CDS_${FEATURE},five_prime_UTR,five_prime_UTR_${FEATURE},three_prime_UTR,three_prime_UTR_${FEATURE}" \
    > ${GENOME_NAME}.trans.${FEATURE}.csv
cat concat.csv \
    | cut -d ',' -f 1,3,5,8,10,13,15,18,20,23,25 \
    >> ${GENOME_NAME}.trans.${FEATURE}.csv
rm concat.csv

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

    runlist stat paralog.yml -s chr.sizes

    # Convert gff3 to runlists
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/feature

    # real	50m49.994s
    # time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
    #     --file Arabidopsis_thaliana.TAIR10.29.gff3.gz \
    #     --size chr.sizes \
    #     --range 2000 \
    #     --clean

    # real	0m41.121s
    time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
        --file ../data/gff3.gz \
        --size ../data/chr.sizes \
        --range 2000
    ```

2. Gene-paralog stats

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha paralog

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha paralog
    ```

3. MITE

    ```bash
    GENOME_NAME=Atha
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

    faops size mite.fa \
        | perl -nla -F"\t" -e '
            next unless defined $F[1];
            $min = defined $min && $min < $F[1] ? $min : $F[1];
            $max = defined $max && $max > $F[1] ? $max : $F[1];
            $n++;
            END{ print qq{seq_number\t$n}; print qq{min_length\t$min}; print qq{max_length\t$max} }' \
        > mite_stat.txt

    find ~/data/alignment/Ensembl/${GENOME_NAME} -type f -name "*.fa" \
        | sort | xargs cat \
        | perl -nl -e '/^>/ or $_ = uc; print' \
        > genome.fa

    perl ~/Scripts/egas/sparsemem_exact.pl -f mite.fa -l 30 -g genome.fa -o mite.replace.tsv
    cat mite.replace.tsv \
        | perl -nla -F"\t" -e ' print for @F' \
        | grep ':' \
        | sort | uniq \
        > mite.position.txt

    wc -l mite.position.txt >> mite_stat.txt
    runlist covers mite.position.txt -o mite.yml

    runlist stat mite.yml -s chr.sizes
    ```

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
    GENOME_NAME=OsatJap
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

    faops size mite.fa \
        | perl -nla -F"\t" -e '
            next unless defined $F[1];
            $min = defined $min && $min < $F[1] ? $min : $F[1];
            $max = defined $max && $max > $F[1] ? $max : $F[1];
            $n++;
            END{ print qq{seq_number\t$n}; print qq{min_length\t$min}; print qq{max_length\t$max} }' \
        > mite_stat.txt

    find ~/data/alignment/Ensembl/${GENOME_NAME} -type f -name "*.fa" \
        | sort | xargs cat \
        | perl -nl -e '/^>/ or $_ = uc; print' \
        > genome.fa

    perl ~/Scripts/egas/sparsemem_exact.pl -f mite.fa -l 25 -g genome.fa -o mite.replace.tsv
    cat mite.replace.tsv \
        | perl -nla -F"\t" -e ' print for @F' \
        | grep ':' \
        | sort | uniq \
        > mite.position.txt

    wc -l mite.position.txt >> mite_stat.txt
    runlist covers mite.position.txt -o mite.yml
    ```

    ```bash
    cd ~/data/alignment/gene-paralog/OsatJap/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh OsatJap mite

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh OsatJap mite
    ```
