# Gene with self-aligning information

## Sources

* [gff3 files](https://github.com/wang-q/withncbi/blob/master/ensembl/ensembl.md#gff3)
* [self-aligning](https://github.com/wang-q/withncbi/blob/master/pop/OPs-selfalign.md#arabidopsis)

## MITE

* http://pmite.hzau.edu.cn/download_mite/
* http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/

## Other repeats

Ensembl gff3 files contain correct descriptions for dust and trf, but repeatmasker's descriptions
is not usable.

So I rerun RepeatMasker on every genomes and get reports from `genome.fa.out`.

### RepeatMasker

Repeat families listed in `genome.fa.tbl`. Families with proportions less than 0.0005 were dropped.

RepeatMasker runned with `-species Viridiplantae`.

* DNA: DNA transposons
* LINE
* LTR
* Low_complexity
* RC: Rolling-circles
* SINE
* Satellite
* Simple_repeat

### Ensembl gff3 repeats

* dust: Low-complexity regions
* trf: Tandem repeats

## Scripts

Same for each species.

### `proc_prepare.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/proc_prepare.sh
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

echo "====> get genome"
find ~/data/alignment/Ensembl/${GENOME_NAME} -type f -name "*.fa" \
    | sort | xargs cat \
    | perl -nl -e '/^>/ or $_ = uc; print' \
    > genome.fa

echo "====> run RepeatMasker"
RepeatMasker genome.fa -species Viridiplantae -xsmall --parallel 8
rm genome.fa.{cat.gz,masked}
rm -fr RM_*

echo "====> paralog stat"
runlist stat paralog.yml -s chr.sizes
mv paralog.yml.csv ../stat

echo "====> Convert gff3 to runlists"
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/feature

# For Atha, 0m41.121s. With --clean, 50m49.994s
time perl ~/Scripts/withncbi/paralog/gff2runlist.pl \
    --file ../data/gff3.gz \
    --size ../data/chr.sizes \
    --range 2000

EOF

```

### `proc_all_gene.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/proc_all_gene.sh
#!/bin/bash

USAGE="Usage: $0 GENOME_NAME FEATURE_FILE"

if [ "$#" -lt 2 ]; then
    echo "$USAGE"
    exit 1
fi

echo "====> parameters"
echo "    " $@

GENOME_NAME=$1
FEATURE_FILE=$2
FEATURE_BASE=`basename "${FEATURE_FILE%.*}"`

cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

runlist stat2 ../feature/all-gene.yml ${FEATURE_FILE} -s ../data/chr.sizes \
    --op intersect --mk --all \
    -o ${GENOME_NAME}.all-gene.${FEATURE_BASE}.csv

EOF

```

### `proc_sep_gene.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/proc_sep_gene.sh
#!/bin/bash

USAGE="Usage: $0 GENOME_NAME FEATURE_FILE"

if [ "$#" -lt 2 ]; then
    echo "$USAGE"
    exit 1
fi

echo "====> parameters"
echo "    " $@

GENOME_NAME=$1
FEATURE_FILE=$2
FEATURE_BASE=`basename "${FEATURE_FILE%.*}"`

cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

echo "====> intersect"
for ftr in gene upstream downstream exon CDS intron five_prime_UTR three_prime_UTR
do
    echo ${ftr}
done \
    | parallel -j 8 --keep-order "
        echo \"====> {} ${FEATURE_BASE}\";
        sleep 1;
        runlist stat2 ../feature/sep-{}.yml ${FEATURE_FILE} -s ../data/chr.sizes \
            --op intersect --mk --all \
            -o stat.sep-{}.${FEATURE_BASE}.csv;
        cat stat.sep-{}.${FEATURE_BASE}.csv \
            | cut -d ',' -f 1,3,5 \
            > stat.sep-{}.${FEATURE_BASE}.csv.tmp;
    "

echo "====> concat gene"
printf "gene_id," > ${GENOME_NAME}.gene.${FEATURE_BASE}.csv
for ftr in gene upstream downstream
do
    printf "${ftr}_length,${ftr}_${FEATURE_BASE},"
done >> ${GENOME_NAME}.gene.${FEATURE_BASE}.csv
echo >> ${GENOME_NAME}.gene.${FEATURE_BASE}.csv

for ftr in gene upstream downstream
do
    cat stat.sep-${ftr}.${FEATURE_BASE}.csv.tmp
done \
    | grep -v "^key" \
    | perl ~/Scripts/withncbi/taxon/merge_csv.pl --concat -f 0 -o stdout \
    >> ${GENOME_NAME}.gene.${FEATURE_BASE}.csv

echo "====> concat trans"
printf "trans_id," > ${GENOME_NAME}.trans.${FEATURE_BASE}.csv
for ftr in exon CDS intron five_prime_UTR three_prime_UTR
do
    printf "${ftr}_length,${ftr}_${FEATURE_BASE},"
done >> ${GENOME_NAME}.trans.${FEATURE_BASE}.csv
echo >> ${GENOME_NAME}.trans.${FEATURE_BASE}.csv

for ftr in exon CDS intron five_prime_UTR three_prime_UTR
do
    cat stat.sep-${ftr}.${FEATURE_BASE}.csv.tmp
done \
    | grep -v "^key" \
    | perl ~/Scripts/withncbi/taxon/merge_csv.pl --concat -f 0 -o stdout \
    >> ${GENOME_NAME}.trans.${FEATURE_BASE}.csv

echo "====> clean"
rm stat.sep-*.${FEATURE_BASE}.csv.tmp
rm stat.sep-*.${FEATURE_BASE}.csv

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

runlist span --op excise -n 10 mite.yml   -o mite.1.yml # remove small spans
runlist span --op fill   -n 10 mite.1.yml -o mite.2.yml # fill small holes
mv mite.2.yml mite.yml
rm mite.1.yml

runlist stat mite.yml -s chr.sizes

echo "==> clean"
rm mite.all.fasta mite.bg.blast mite.bg.fasta mite.filter.fa mite.replace.tsv

mv mite.yml.csv ../stat
mv mite_stat.txt ../stat

EOF

```

### `proc_repeat.sh`

```bash

cat <<'EOF' > ~/data/alignment/gene-paralog/proc_repeat.sh
#!/bin/bash

USAGE="Usage: $0 GENOME_NAME"

if [ "$#" -lt 1 ]; then
    echo "$USAGE"
    exit 1
fi

echo "====> parameters"
echo "    " $@

GENOME_NAME=$1
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat

find . -type f -name "*.yml" -or -name "*.csv" | parallel rm

echo "==> gff types"
gzip -d -c ../data/gff3.gz \
    | perl -nla -e '/^#/ and next; print $F[2]' \
    | sort | uniq -c \
    > gff.type.txt

gzip -d -c ../data/gff3.gz \
    | perl -nla -e '
        /^#/ and next;
        $F[2] eq q{repeat_region} or next;
        $F[8] =~ /description\=(\w+)/i or next;
        print qq{$F[2]\t$F[1]\t$1};
        ' \
    | sort | uniq -c \
    > gff.repeat.txt

echo "==> rmout families"
cat ../data/genome.fa.out \
    | perl -nla -e '/^\s*\d+/ or next; print $F[10]' \
    | sort | uniq -c \
    > rmout.family.txt

echo "==> rmout results"
perl ~/Scripts/withncbi/paralog/rmout2runlist.pl \
    --file ../data/genome.fa.out \
    --size ../data/chr.sizes

runlist split ../feature/all-repeat.yml -o .
runlist merge *.yml -o all-repeat.yml
find . -type f -name "*.yml" -not -name "all-*" | parallel rm

runlist span --op excise -n 10 --mk all-repeat.yml -o all-repeat.1.yml   # remove small spans
runlist span --op fill   -n 10 --mk all-repeat.1.yml -o all-repeat.2.yml # fill small holes
runlist split all-repeat.2.yml -o .
mv all-repeat.2.yml all-repeat.yml
rm all-repeat.1.yml

echo "==> basic stat"
runlist stat all-repeat.yml -s ../data/chr.sizes --mk --all

echo "==> find repeat families large enough"
cat all-repeat.yml.csv \
    | perl -nla -F"," -e '
        next if $F[3] =~ /^c/;
        print $F[0] if $F[3] > 0.0005
    ' \
    > repeat.family.txt

echo "==> clean"
mv all-repeat.yml.csv ../stat

EOF

```

## Atha

1. Prepare

    ```bash
    GENOME_NAME=Atha

    echo "====> create directories"
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

    echo "====> copy or download needed files here"
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    cp ~/data/ensembl82/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.29.gff3.gz gff3.gz
    cp ~/data/alignment/self/arabidopsis/Genomes/Atha/chr.sizes chr.sizes
    cp ~/data/alignment/self/arabidopsis/Results/Atha/Atha.chr.runlist.yml paralog.yml
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/03_arabidopsis_mite_seq.fa -O mite.fa
    ```

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/data

    bash ~/data/alignment/gene-paralog/proc_prepare.sh Atha
    ```

2. Gene-paralog stats

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha ../data/paralog.yml

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha ../data/paralog.yml
    ```

3. MITE

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/data

    bash ~/data/alignment/gene-paralog/proc_mite.sh Atha
    ```

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat

    bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha ../data/mite.yml

    bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha ../data/mite.yml
    ```

4. Other Repeats

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/repeat

    bash ~/data/alignment/gene-paralog/proc_repeat.sh Atha
    ```

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat

    cat ../repeat/repeat.family.txt \
        | parallel -j 8 --keep-order "
            bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha ../repeat/{}.yml
        "

    cat ../repeat/repeat.family.txt \
        | parallel -j 1 --keep-order "
            bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha ../repeat/{}.yml
        "
    ```

## Plants aligned with full chromosomes

1. Prepare

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        echo "====> create directories"
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

        echo "====> copy or download needed files here"
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
        cp ~/data/alignment/self/plants/Genomes/${GENOME_NAME}/chr.sizes chr.sizes
        cp ~/data/alignment/self/plants/Results/${GENOME_NAME}/${GENOME_NAME}.chr.runlist.yml paralog.yml
    done

    # http://stackoverflow.com/questions/1494178/how-to-define-hash-tables-in-bash
    # OSX has bash 3. So no easy hashmaps. Do it manually.
    cd ~/data/alignment/gene-paralog

    # OsatJap
    cp ~/data/ensembl82/gff3/oryza_sativa/Oryza_sativa.IRGSP-1.0.29.gff3.gz \
        ~/data/alignment/gene-paralog/OsatJap/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/26_nipponbare_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/OsatJap/data/mite.fa

    # Alyr
    cp ~/data/ensembl82/gff3/arabidopsis_lyrata/Arabidopsis_lyrata.v.1.0.29.gff3.gz \
        ~/data/alignment/gene-paralog/Alyr/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/02_lyrata_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Alyr/data/mite.fa

    # Sbic
    # Ensemblgenomes 82 didn't provide a full gff3
    gt merge -gzip -force \
        -o ~/data/alignment/gene-paralog/Sbic/data/gff3.gz \
        ~/data/ensembl82/gff3/sorghum_bicolor/Sorghum_bicolor.Sorbi1.29.chromosome.*.gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/28_sorghum_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Sbic/data/mite.fa

    ```

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
        bash ~/data/alignment/gene-paralog/proc_prepare.sh ${GENOME_NAME}
    done
    ```

2. Gene-paralog stats

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        bash ~/data/alignment/gene-paralog/proc_all_gene.sh ${GENOME_NAME} ../data/paralog.yml
        bash ~/data/alignment/gene-paralog/proc_sep_gene.sh ${GENOME_NAME} ../data/paralog.yml
    done
    ```

3. MITE

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
        bash ~/data/alignment/gene-paralog/proc_mite.sh ${GENOME_NAME}
    done
    ```

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        bash ~/data/alignment/gene-paralog/proc_all_gene.sh ${GENOME_NAME} ../data/mite.yml
        bash ~/data/alignment/gene-paralog/proc_sep_gene.sh ${GENOME_NAME} ../data/mite.yml
    done
    ```

4. Other Repeats

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat
        bash ~/data/alignment/gene-paralog/proc_repeat.sh ${GENOME_NAME}
    done
    ```

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
	    cat ../repeat/repeat.family.txt \
	        | parallel -j 8 --keep-order "
	            bash ~/data/alignment/gene-paralog/proc_all_gene.sh ${GENOME_NAME} ../repeat/{}.yml
	        "

	    cat ../repeat/repeat.family.txt \
	        | parallel -j 1 --keep-order "
	            bash ~/data/alignment/gene-paralog/proc_sep_gene.sh ${GENOME_NAME} ../repeat/{}.yml
	        "
    done
    ```
