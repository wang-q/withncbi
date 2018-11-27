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

# Atha

Full processing time is about 1 hour.

1. [Data](https://github.com/wang-q/withncbi/blob/master/paralog/OPs-selfalign.md#arabidopsis)

2. Prepare

    ```bash
    GENOME_NAME=Atha

    echo "====> create directories"
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
    mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/yml

    echo "====> copy or download needed files here"
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data
    cp ~/data/alignment/self/plants/Genomes/${GENOME_NAME}/chr.sizes chr.sizes
    cp ~/data/alignment/self/plants/Results/${GENOME_NAME}/${GENOME_NAME}.chr.runlist.yml paralog.yml

    cp ~/data/ensembl82/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.29.gff3.gz gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/03_arabidopsis_mite_seq.fa -O mite.fa
    ```

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/data

    # 0m44.430s
    time bash ~/data/alignment/gene-paralog/proc_prepare.sh Atha
    # 0m46.052s
    time bash ~/data/alignment/gene-paralog/proc_repeat.sh Atha
    # 0m43.666s
    time bash ~/data/alignment/gene-paralog/proc_mite.sh Atha
    ```

3. Paralog-repeats stats

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat
    # 0m28.356s
    time bash ~/data/alignment/gene-paralog/proc_paralog.sh Atha
    ```

4. Gene-paralog stats

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat
    # 0m6.888s
    time bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha ../yml/paralog.yml
    # 0m7.283s
    time bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha ../yml/paralog_adjacent.yml

    # E5-2690 v3
    # real    8m45.668s
    # user    57m58.984s
    # sys     0m1.808s
    # i7-6700k
    # real	15m18.045s
    # user	104m21.728s
    # sys	0m13.363s
    time bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha ../yml/paralog.yml

    # real	0m35.566s
    # user	1m12.613s
    # sys	0m7.601s
    time bash ~/data/alignment/gene-paralog/proc_sep_gene_jrunlist.sh Atha ../yml/paralog.yml

    bash ~/data/alignment/gene-paralog/proc_sep_gene_jrunlist.sh Atha ../yml/paralog_adjacent.yml
    ```

5. Gene-repeats stats

    ```bash
    cd ~/data/alignment/gene-paralog/Atha/stat

    cat ../yml/repeat.family.txt \
        | parallel -j 8 --keep-order "
            bash ~/data/alignment/gene-paralog/proc_all_gene.sh Atha ../yml/{}.yml
        "

    # 12 hours?
    # time \
    # cat ../yml/repeat.family.txt \
    #     | parallel -j 1 --keep-order "
    #         bash ~/data/alignment/gene-paralog/proc_sep_gene.sh Atha ../yml/{}.yml
    #     "    

    # real	10m34.669s
    # user	17m27.433s
    # sys	3m20.765s
    time \
        cat ../yml/repeat.family.txt \
            | parallel -j 1 --keep-order "
                bash ~/data/alignment/gene-paralog/proc_sep_gene_jrunlist.sh Atha ../yml/{}.yml
            "
    ```

6. Pack up

    ```bash
    cd ~/data/alignment/gene-paralog
    find Atha -type f -not -path "*/data/*" -print | zip Atha.zip -9 -@
    ```

# Plants aligned with full chromosomes

1. [Data](https://github.com/wang-q/withncbi/blob/master/paralog/OPs-selfalign.md#full-chromosomes)

    * OsatJap
    * Bdis
    * Alyr
    * Sbic

2. Prepare

    ```bash
    for GENOME_NAME in OsatJap Bdis Alyr Sbic
    do
        echo "====> create directories"
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/data
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/feature
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/yml

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

    # Bdis
    cp ~/data/ensembl82/gff3/brachypodium_distachyon/Brachypodium_distachyon.v1.0.29.gff3.gz \
        ~/data/alignment/gene-paralog/Bdis/data/gff3.gz
    wget http://pmite.hzau.edu.cn/MITE/MITE-SEQ-V2/25_brachypodium_mite_seq.fa \
        -O ~/data/alignment/gene-paralog/Bdis/data/mite.fa

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

        bash ~/data/alignment/gene-paralog/proc_repeat.sh ${GENOME_NAME}
        bash ~/data/alignment/gene-paralog/proc_mite.sh ${GENOME_NAME}
    done
    ```

3. Paralog-repeats stats

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        bash ~/data/alignment/gene-paralog/proc_paralog.sh ${GENOME_NAME}
    done
    ```

4. Gene-paralog stats

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
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
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
        cat ../yml/repeat.family.txt \
            | parallel -j 8 --keep-order "
                bash ~/data/alignment/gene-paralog/proc_all_gene.sh ${GENOME_NAME} ../yml/{}.yml
            "

        cat ../yml/repeat.family.txt \
            | parallel -j 1 --keep-order "
                bash ~/data/alignment/gene-paralog/proc_sep_gene.sh ${GENOME_NAME} ../yml/{}.yml
            "
    done
    ```

6. Pack up

    ```bash
    for GENOME_NAME in OsatJap Alyr Sbic
    do
        cd ~/data/alignment/gene-paralog
        find ${GENOME_NAME} -type f -not -path "*/data/*" -print | zip ${GENOME_NAME}.zip -9 -@
    done
    ```

# Plants aligned with partitioned chromosomes

1. [Data](https://github.com/wang-q/withncbi/blob/master/paralog/OPs-selfalign.md#partitioned-chromosomes)

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

