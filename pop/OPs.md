# Operating steps for each group.

Less detailed than Trichoderma in [README.md](README.md), but include examples
for genomes out of WGS, which usually in better assembling levels.

## *Saccharomyces* WGS

1. Create `pop/saccharomyces.tsv` manually.

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Saccharomyces&order=organism
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid4930[Organism:exp]

    ```bash
    export GENUS_ID=4930
    export GENUS=saccharomyces
    mkdir -p ~/data/alignment/$GENUS
    cd ~/data/alignment/$GENUS

    ...

    unset GENUS_ID
    unset GENUS
    ```

    Remove all strains of Saccharomyces cerevisiae.

    ```bash
    cat saccharomyces.tsv | grep Scer_ > scer_.tsv
    cat saccharomyces.tsv | grep -v Scer_ > sacc_.tsv
    ```

    Edit them to fix names and comment out bad strains.

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/saccharomyces
    cd ~/data/alignment/saccharomyces

    perl ~/Scripts/withncbi/util/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/saccharomyces.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -x 6 -s 3 -c -i WGS/saccharomyces.url.txt

    find WGS -name "*.gz" | xargs gzip -t

    # rsync remote files
    # rsync --progress -av wangq@139.162.23.84:/home/wangq/data/alignment/saccharomyces/ ~/data/alignment/saccharomyces
    ```

3. Download *Saccharomyces cerevisiae* S288c

    ```bash
    cd ~/data/alignment/saccharomyces/WGS

    # Omit chrMt
    perl ~/Scripts/withncbi/util/assemble_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt \
        -name Scer_S288c \
        --nuclear \
        > Scer_S288c.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/util/batch_get_seq.pl \
        -f Scer_S288c.seq.csv  \
        -r -p
    ```

4. `gen_pop_conf.pl`

    ```bash
    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/saccharomyces/WGS/saccharomyces.data.yml \
        -o ~/Scripts/withncbi/pop/saccharomyces_test.yml \
        -d ~/data/alignment/saccharomyces/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=saccharomyces \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/saccharomyces' \
        --opt rm_species=Fungi \
        --downloaded 'name=Scer_S288c;taxon=559292;sciname=Saccharomyces cerevisiae S288c' \
        --plan 'name=four_way;t=Scer_S288c;qs=Sbou_ATCC_MYA_796,Spar_NRRL_Y_17217,Spas_CBS_1483'
    ```

5. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/saccharomyces_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_four_way.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

## *Candida* WGS

1. Create `pop/candida.tsv` manually.

    Check NCBI pages

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Candida*&order=organism
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid1535326[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid1535326[Organism:exp]

    The wgs page contains a lot of strains from other genus.

    ```bash
    export GENUS_ID=1535326
    export GENUS=candida
    mkdir -p ~/data/alignment/$GENUS
    cd ~/data/alignment/$GENUS

    ...
    # Edit raw2.tsv, remove lines containing CANDIDATUS or CANDIDATE DIVISION
    cat raw2.tsv | grep -v 'CANDIDATUS' | grep -v 'CANDIDATE DIVISION' > tmp.tsv
    mv tmp.tsv raw2.tsv
    ...

    unset GENUS_ID
    unset GENUS
    ```


2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/candida
    cd ~/data/alignment/candida

    perl ~/Scripts/withncbi/util/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/candida.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -x 6 -s 3 -c -i WGS/candida.url.txt

    find WGS -name "*.gz" | xargs gzip -t

    # rsync remote files
    # rsync --progress -av wangq@139.162.23.84:/home/wangq/data/alignment/candida/ ~/data/alignment/candida
    ```

3. Download *Candida dubliniensis* CD36 and *Candida orthopsilosis* Co 90-125

    ```bash
    cd ~/data/alignment/candida/WGS

    perl ~/Scripts/withncbi/util/assemble_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000026945.1.assembly.txt \
        -name Cdub_CD36 \
        > Cdub_CD36.seq.csv

    perl ~/Scripts/withncbi/util/assemble_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000315875.1.assembly.txt \
        -name Corh_Co_90_125 \
        > Corh_Co_90_125.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/util/batch_get_seq.pl \
        -f Cdub_CD36.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/util/batch_get_seq.pl \
        -f Corh_Co_90_125.seq.csv  \
        -r -p
    ```

4. `gen_pop_conf.pl`

    Pay attentions to --downloaded orders. The first one will be the default target.

    ```bash
    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/candida/WGS/candida.data.yml \
        -o ~/Scripts/withncbi/pop/candida_test.yml \
        -d ~/data/alignment/candida/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=candida \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/candida' \
        --opt rm_species=Fungi \
        --downloaded 'name=Cdub_CD36;taxon=573826;sciname=Candida dubliniensis CD36' \
        --downloaded 'name=Corh_Co_90_125;taxon=1136231;sciname=Candida orthopsilosis Co 90-125' \
        --plan 'name=four_way;t=Cdub_CD36;qs=Corh_Co_90_125,Calb_WO_1,Ctro_MYA_3404' \
        --plan 'name=four_way_2;t=Corh_Co_90_125;qs=Cdub_CD36,Calb_WO_1,Ctro_MYA_3404'
    ```

5. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/candida_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_four_way.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_four_way_2.sh

    sh 3_pair_cmd.sh # Only do this when target switched
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

## *Fusarium* WGS

1. Create `pop/fusarium.tsv` manually.

    Check NCBI pages

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=fusarium*&order=organism
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid5506[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid5506[Organism:exp]

    This genus contains about 40 strains, use command lines listed in README.md.

    ```bash
    export GENUS_ID=5506
    export GENUS=fusarium
    mkdir -p ~/data/alignment/$GENUS
    cd ~/data/alignment/$GENUS

    ...

    unset GENUS_ID
    unset GENUS
    ```

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/fusarium
    cd ~/data/alignment/fusarium

    perl ~/Scripts/withncbi/util/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/fusarium.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -x 6 -s 3 -c -i WGS/fusarium.url.txt

    find WGS -name "*.gz" | xargs gzip -t

    # rsync remote files
    # rsync --progress -av wangq@139.162.23.84:/home/wangq/data/alignment/fusarium/ ~/data/alignment/fusarium
    ```

3. Download Fusarium graminearum* PH-1, *Fusarium pseudograminearum* CS3270
  and *Fusarium verticillioides* 7600.

    ```bash
    cd ~/data/alignment/candida/WGS

    perl ~/Scripts/withncbi/util/assemble_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000240135.3.assembly.txt \
        -name Fgra_PH_1 \
        --genbank \
        --nuclear \
        > Fgra_PH_1.seq.csv

    perl ~/Scripts/withncbi/util/assemble_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000974265.1.assembly.txt \
        -name Fpse_CS3270 \
        --genbank \
        --nuclear \
        > Fpse_CS3270.seq.csv

    perl ~/Scripts/withncbi/util/assemble_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000149555.1.assembly.txt \
        -name Fver_7600 \
        --genbank \
        --nuclear \
        > Fver_7600.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/util/batch_get_seq.pl \
        -f Fgra_PH_1.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/util/batch_get_seq.pl \
        -f Fpse_CS3270.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/util/batch_get_seq.pl \
        -f Fver_7600.seq.csv  \
        -r -p
    ```

3. `gen_pop_conf.pl`

    Pay attentions to --downloaded orders. The first one will be the default target.

    ```bash
    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/fusarium/WGS/fusarium.data.yml \
        -o ~/Scripts/withncbi/pop/fusarium_test.yml \
        -d ~/data/alignment/fusarium/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=fusarium \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/fusarium' \
        --opt rm_species=Fungi \
        --downloaded 'name=Cdub_CD36;taxon=573826;sciname=Candida dubliniensis CD36' \
        --downloaded 'name=Corh_Co_90_125;taxon=1136231;sciname=Candida orthopsilosis Co 90-125' \
        --plan 'name=four_way;t=Cdub_CD36;qs=Corh_Co_90_125,Calb_WO_1,Ctro_MYA_3404' \
        --plan 'name=four_way_2;t=Corh_Co_90_125;qs=Cdub_CD36,Calb_WO_1,Ctro_MYA_3404'
    ```

4. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/fusarium_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_XXX.sh

    # sh 3_pair_cmd.sh # Only do this when target switched, e.g. four_way_2
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```
