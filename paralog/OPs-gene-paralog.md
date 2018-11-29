# Gene-paralog-repeats

[TOC levels=1-3]: # " "
- [Gene-paralog-repeats](#gene-paralog-repeats)
- [Sources](#sources)
- [TODO](#todo)
- [Stats](#stats)
- [Repeats](#repeats)
    - [MITE](#mite)
    - [Other repeats](#other-repeats)
- [Scripts](#scripts)
- [Create directories](#create-directories)
- [Symlink ensembl gff files](#symlink-ensembl-gff-files)
- [`chr.sizes` and `paralog.yml`](#chrsizes-and-paralogyml)
- [Download MITE](#download-mite)
- [Atha](#atha)
- [Plants aligned with full chromosomes](#plants-aligned-with-full-chromosomes)
- [Plants aligned with partitioned chromosomes](#plants-aligned-with-partitioned-chromosomes)
- [Plants with annotations from JGI](#plants-with-annotations-from-jgi)


# Sources

* Annotations from
  [Ensembl gff3 files](https://github.com/wang-q/withncbi/blob/master/ensembl/README.md#eg-gff3)
* Annotations from JGI PhytozomeV12.1
    * `PhytozomeV12_unrestricted/Athaliana`
    * `PhytozomeV12_unrestricted/Osativa`
* Paralogs from
  [self-aligning](https://github.com/wang-q/withncbi/blob/master/paralog/OPs-selfalign.md)

# TODO

* Remove full transposons (Retro, DNA and RC transposons) from `paralog.yml`.
* Subfamilies.

# Stats

* Coverages on chromosomes of all feature types.
* Paralogs and adjacent regions intersect with all repeat families.
    * paralog
    * paralog_adjacent: paralog + 2000 bp
    * paralog_gene: intersections between paralogs and gene + 2000 bp
* Genes, upstreams, downstreams intersect with paralogs and all repeat families.
    * Up/down-streams are 2000 bp.
* Exons, introns, CDSs, five_prime_UTRs and three_prime_UTRs intersect with paralogs and all repeat
  families.

# Repeats

## MITE

* [Plant MITE database](http://pmite.hzau.edu.cn/download_mite/)
* http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/

## Other repeats

Ensembl gff3 files contain correct descriptions for dust and trf, but repeatmasker's descriptions
are not usable.

So I rerun RepeatMasker on every genomes and get reports from `genome.fa.out`.

1. RepeatMasker

    Run RepeatMasker with `-species Viridiplantae`.

    Repeat families listed in `genome.fa.tbl`. Families with proportions less than **0.0005** were
    dropped.

    * DNA: DNA transposons
    * LINE
    * LTR
    * Low_complexity
    * RC: Rolling-circles
    * SINE
    * Satellite
    * Simple_repeat

2. Ensembl gff3 repeats

    * dust: Low-complexity regions
    * trf: Tandem repeats

# Scripts

Same for each species.

* `proc_prepare.sh`

    Genome, RepeatMasker, dustmasker and gff3.

* `proc_repeat.sh`

    Repeats from RepeatMasker and gff3. Create `repeat.family.txt`.

* `proc_mite.sh`

    MITE. Append to `repeat.family.txt`.

* `proc_paralog.sh`

* `proc_all_gene.sh`

* `proc_sep_gene.sh`

# Create directories

```bash
for GENOME_NAME in Atha Alyr OsatJap Sbic; do
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/yml
done

```

# Symlink ensembl gff files

```bash
for GENOME_NAME in Atha Alyr OsatJap Sbic; do
    echo "==> ${GENOME_NAME}"
    
    pushd ~/data/alignment/gene-paralog/${GENOME_NAME}/data > /dev/null
    ln -s ~/data/alignment/Ensembl/${GENOME_NAME}/chr.gff chr.gff
    popd > /dev/null
    
    echo
done

```

```bash
# Sbic
## Ensemblgenomes 82 didn't provide a full gff3
#gt merge -gzip -force \
#    -o ~/data/alignment/gene-paralog/Sbic/data/gff3.gz \
#    ~/data/ensembl82/gff3/sorghum_bicolor/Sorghum_bicolor.Sorbi1.29.chromosome.*.gff3.gz

```

# `chr.sizes` and `paralog.yml`

```bash
for GENOME_NAME in Atha Alyr OsatJap Sbic; do
    echo "==> ${GENOME_NAME}"
    
    pushd ~/data/alignment/gene-paralog/${GENOME_NAME}/data > /dev/null
    cp ~/data/alignment/self/plants/Processing/${GENOME_NAME}/chr.sizes chr.sizes
    cp ~/data/alignment/self/plants/Results/${GENOME_NAME}/${GENOME_NAME}.cover.yml paralog.yml
    popd > /dev/null
    
    echo
done

```

# Download MITE

```bash

# http://stackoverflow.com/questions/1494178/how-to-define-hash-tables-in-bash
# https://stackoverflow.com/questions/6047648/bash-4-associative-arrays-error-declare-a-invalid-option

find ~/data/alignment/gene-paralog/ -name "mite.fa" | xargs rm

ARRAY=(
    'Atha::http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/03_arabidopsis_mite_seq.fa'
    'Alyr::http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/02_lyrata_mite_seq.fa'
    'OsatJap::http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/26_nipponbare_mite_seq.fa'
    'Sbic::http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/28_sorghum_mite_seq.fa'
)

for item in "${ARRAY[@]}" ; do
    GENOME_NAME="${item%%::*}"
    MITE_URL="${item##*::}"
    FILENAME=$(MITE_URL="${MITE_URL}" perl -e '$ENV{MITE_URL} =~ m{\/(\w+\.fa)$}; print $1')
    
    echo "==> ${FILENAME}"

    pushd ~/data/alignment/gene-paralog/${GENOME_NAME}/data > /dev/null
    wget -N "${MITE_URL}"
    ln -s ${FILENAME} mite.fa
    popd > /dev/null
    
    echo
done

```

# Atha

Full processing time is about 1 hour.

* [Data](OPs-selfalign.md#arabidopsis)

* Prepare

```bash
cd ~/data/alignment/gene-paralog/Atha/data

bash ~/Scripts/withncbi/paralog/proc_prepare.sh Atha
bash ~/Scripts/withncbi/paralog/proc_repeat.sh Atha
bash ~/Scripts/withncbi/paralog/proc_mite.sh Atha
```

* Paralog-repeats stats

```bash
cd ~/data/alignment/gene-paralog/Atha/stat
bash ~/Scripts/withncbi/paralog/proc_paralog.sh Atha
```

* Gene-paralog stats

```bash
cd ~/data/alignment/gene-paralog/Atha/stat
bash ~/Scripts/withncbi/paralog/proc_all_gene.sh Atha ../yml/paralog.yml
bash ~/Scripts/withncbi/paralog/proc_all_gene.sh Atha ../yml/paralog_adjacent.yml

time bash ~/Scripts/withncbi/paralog/proc_sep_gene.sh Atha ../yml/paralog.yml
#real    1m57.614s
#user    1m27.512s
#sys     2m3.066s

time bash ~/Scripts/withncbi/paralog/proc_sep_gene.sh Atha ../yml/paralog_adjacent.yml
#real    2m31.327s
#user    1m34.904s
#sys     3m13.446s

```

* Gene-repeats stats

```bash
cd ~/data/alignment/gene-paralog/Atha/stat

time cat ../yml/repeat.family.txt |
    parallel -j 4 --keep-order '
        echo "==> {}"
        bash ~/Scripts/withncbi/paralog/proc_all_gene.sh Atha ../yml/{}.yml
        echo
    '
#real    0m19.571s
#user    1m2.364s
#sys     0m0.925s

time cat ../yml/repeat.family.txt |
    parallel -j 1 --keep-order '
        echo "==> {}"
        bash ~/Scripts/withncbi/paralog/proc_sep_gene.sh Atha ../yml/{}.yml
        echo
    '
#real    18m21.261s
#user    17m41.651s
#sys     33m40.853s

```

* Pack up

```bash
cd ~/data/alignment/gene-paralog
find Atha -type f -not -path "*/data/*" -print | zip Atha.zip -9 -@
```

# Plants aligned with full chromosomes

* [Data](OPs-selfalign.md#plants-full-chromosomes)

    * Alyr
    * OsatJap
    * Sbic

```bash
# Prepare
for GENOME_NAME in Alyr OsatJap Sbic; do
    echo "==> ${GENOME_NAME}"
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

    bash ~/Scripts/withncbi/paralog/proc_prepare.sh ${GENOME_NAME}

    bash ~/Scripts/withncbi/paralog/proc_repeat.sh ${GENOME_NAME}
    bash ~/Scripts/withncbi/paralog/proc_mite.sh ${GENOME_NAME}
    echo
done

# Paralog-repeats stats
for GENOME_NAME in Alyr OsatJap Sbic; do
    echo "==> ${GENOME_NAME}"
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
    bash ~/Scripts/withncbi/paralog/proc_paralog.sh ${GENOME_NAME}
    echo
done

# Gene-paralog stats
for GENOME_NAME in Alyr OsatJap Sbic; do
    echo "==> ${GENOME_NAME}"
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

    bash ~/Scripts/withncbi/paralog/proc_all_gene.sh ${GENOME_NAME} ../yml/paralog.yml
    bash ~/Scripts/withncbi/paralog/proc_all_gene.sh ${GENOME_NAME} ../yml/paralog_adjacent.yml

    bash ~/Scripts/withncbi/paralog/proc_sep_gene.sh ${GENOME_NAME} ../yml/paralog.yml
    bash ~/Scripts/withncbi/paralog/proc_sep_gene.sh ${GENOME_NAME} ../yml/paralog_adjacent.yml
    echo
done

# Gene-repeats stats
for GENOME_NAME in Alyr OsatJap Sbic; do
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
    cat ../yml/repeat.family.txt |
        parallel -j 4 --keep-order "
            echo '==> {}'
            bash ~/Scripts/withncbi/paralog/proc_all_gene.sh ${GENOME_NAME} ../yml/{}.yml
            echo
        "

    cat ../yml/repeat.family.txt |
        parallel -j 1 --keep-order "
            echo '==> {}'
            bash ~/Scripts/withncbi/paralog/proc_sep_gene.sh ${GENOME_NAME} ../yml/{}.yml
            echo
        "
done

# Pack up
for GENOME_NAME in Alyr OsatJap Sbic; do
    echo "==> ${GENOME_NAME}"
    cd ~/data/alignment/gene-paralog
    find ${GENOME_NAME} -type f -not -path "*/data/*" -print | zip ${GENOME_NAME}.zip -9 -@
    echo
done

```

# Plants aligned with partitioned chromosomes

1. [Data](OPs-selfalign.md#plants-partitioned-chromosomes)

    * Mtru
    * Gmax
    * Brap
    * Vvin
    * Slyc
    * Stub

    * Bole (no mite)

2. Prepare

    ```bash
    for GENOME_NAME in Mtru Gmax Brap Vvin Slyc Stub
    do
        echo "====> create directories"
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/yml

        echo "====> copy or download needed files here"
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
        cp ~/data/alignment/self/plants_parted/Genomes/${GENOME_NAME}/chr.sizes chr.sizes
        cp ~/data/alignment/self/plants_parted/Results/${GENOME_NAME}/${GENOME_NAME}.chr.runlist.yml paralog.yml
    done

    cd ~/data/alignment/gene-paralog

    # Mtru
    cp ~/data/ensembl82/gff3/medicago_truncatula/Medicago_truncatula.GCA_000219495.2.29.gff3.gz \
        ~/data/alignment/gene-paralog/Mtru/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/20_medicago_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Mtru/data/mite.fa

    # Gmax
    cp ~/data/ensembl82/gff3/glycine_max/Glycine_max.V1.0.29.gff3.gz \
        ~/data/alignment/gene-paralog/Gmax/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/18_soybean_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Gmax/data/mite.fa

    # Brap
    cp ~/data/ensembl82/gff3/brassica_rapa/Brassica_rapa.IVFCAASv1.29.gff3.gz \
        ~/data/alignment/gene-paralog/Brap/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/04_brassica_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Brap/data/mite.fa

    # Vvin
    cp ~/data/ensembl82/gff3/vitis_vinifera/Vitis_vinifera.IGGP_12x.29.gff3.gz \
        ~/data/alignment/gene-paralog/Vvin/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/39_grape_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Vvin/data/mite.fa

    # Slyc
    cp ~/data/ensembl82/gff3/solanum_lycopersicum/Solanum_lycopersicum.GCA_000188115.2.29.gff3.gz \
        ~/data/alignment/gene-paralog/Slyc/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/37_tomato_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Slyc/data/mite.fa

    # Stub
    cp ~/data/ensembl82/gff3/solanum_tuberosum/Solanum_tuberosum.3.0.29.gff3.gz \
        ~/data/alignment/gene-paralog/Stub/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/38_potato_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Stub/data/mite.fa

    # Bdis
    cp ~/data/ensembl82/gff3/brachypodium_distachyon/Brachypodium_distachyon.v1.0.29.gff3.gz \
        ~/data/alignment/gene-paralog/Bdis/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/25_brachypodium_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Bdis/data/mite.fa

    # # Sita
    # cp ~/data/ensembl82/gff3/setaria_italica/Setaria_italica.JGIv2.0.29.gff3.gz \
    #     ~/data/alignment/gene-paralog/Sita/data/gff3.gz
    # wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/27_foxtail_mite_seq.fa \
    #     -O ~/data/alignment/gene-paralog/Sita/data/mite.fa
    ```

    ```bash
    for GENOME_NAME in Mtru Gmax Brap Vvin Slyc Stub Sita
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

        bash ~/data/alignment/gene-paralog/proc_prepare.sh ${GENOME_NAME}

        bash ~/data/alignment/gene-paralog/proc_repeat.sh ${GENOME_NAME}
        bash ~/data/alignment/gene-paralog/proc_mite.sh ${GENOME_NAME}
    done
    ```

# Plants with annotations from JGI

1. [Data](https://github.com/wang-q/withncbi/blob/master/paralog/OPs-selfalign.md#full-chromosomes)

    * AthaJGI
    * OsatJapJGI

2. Prepare

    ```bash
    for GENOME_NAME in Atha OsatJap
    do
        echo "====> create directories"
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}JGI/data
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}JGI/feature
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}JGI/repeat
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}JGI/stat
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}JGI/yml

        echo "====> copy or download needed files here"
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}JGI/data
        cp ~/data/alignment/gene-paralog/${GENOME_NAME}/data/* .
    done

    cd ~/data/alignment/gene-paralog

    # Atha
    cp -f ~/data/PhytozomeV11/Athaliana/annotation/Athaliana_167_TAIR10.gene_exons.gff3.gz \
        ~/data/alignment/gene-paralog/AthaJGI/data/gff3.gz

    # OsatJap
    cp ~/data/PhytozomeV11/Osativa/annotation/Osativa_323_v7.0.gene_exons.gff3.gz \
        ~/data/alignment/gene-paralog/OsatJapJGI/data/gff3.gz
    ```

    ```bash
    for GENOME_NAME in AthaJGI OsatJapJGI
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
        perl ~/Scripts/withncbi/util/gff2runlist.pl \
            --file ../data/gff3.gz \
            --size ../data/chr.sizes \
            --range 2000 --remove

        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

        #bash ~/data/alignment/gene-paralog/proc_prepare.sh ${GENOME_NAME}

        bash ~/data/alignment/gene-paralog/proc_repeat.sh ${GENOME_NAME}
        bash ~/data/alignment/gene-paralog/proc_mite.sh ${GENOME_NAME}
    done
    ```

3. Paralog-repeats stats

    ```bash
    for GENOME_NAME in AthaJGI OsatJapJGI
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        bash ~/data/alignment/gene-paralog/proc_paralog.sh ${GENOME_NAME}
    done
    ```

4. Gene-paralog stats

    ```bash
    for GENOME_NAME in AthaJGI OsatJapJGI
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

        bash ~/data/alignment/gene-paralog/proc_all_gene.sh ${GENOME_NAME} ../yml/paralog.yml
        bash ~/data/alignment/gene-paralog/proc_all_gene.sh ${GENOME_NAME} ../yml/paralog_adjacent.yml

        bash ~/data/alignment/gene-paralog/proc_sep_gene.sh ${GENOME_NAME} ../yml/paralog.yml
        bash ~/data/alignment/gene-paralog/proc_sep_gene.sh ${GENOME_NAME} ../yml/paralog_adjacent.yml
    done
    ```

5. Gene-repeats stats

    ```bash
    for GENOME_NAME in AthaJGI OsatJapJGI
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        cat ../yml/repeat.family.txt \
            | parallel -j 8 --keep-order "
                bash ~/data/alignment/gene-paralog/proc_all_gene.sh ${GENOME_NAME} ../yml/{}.yml
            "

        cat ../yml/repeat.family.txt \
            | parallel -j 1 --keep-order "
                bash ~/data/alignment/gene-paralog/proc_sep_gene_jrunlist.sh ${GENOME_NAME} ../yml/{}.yml
            "
    done
    ```

6. Pack up

    ```bash
    for GENOME_NAME in AthaJGI OsatJapJGI
    do
        cd ~/data/alignment/gene-paralog
        find ${GENOME_NAME} -type f -not -path "*/data/*" -print | zip ${GENOME_NAME}.zip -9 -@
    done
    ```

