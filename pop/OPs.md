# Operating steps for each group.

## Saccharomyces WGS

1. Create pop/saccharomyces.tsv manually.

    http://www.ncbi.nlm.nih.gov/assembly?term=Saccharomyces%20kudriavzevii

2. Create working direcotry and download WGS sequences.

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

3. Download Saccharomyces cerevisiae S288c

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

3. 'gen_pop_conf.pl`

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
        --plan 'name=four_way;t=Scer_S288c;qs=Sbou_ATCC_MYA_796,Spar_NRRL_Y_17217,Spas_CBS_1483' \
        --plan 'name=17_way;t=Scer_S288c;qs=Sarb_H_6,Sbay_623_6C,Sbou_17,Sbou_ATCC_MYA_796,Sbou_EDRL,ScerSkud_VIN7,Skud_IFO_1802,Skud_ZP591,Smik_IFO_1815_1,Spar_NRRL_Y_17217,Spas_CBS_1483,Spas_CBS_1513,Spas_CCY48_91,Spas_Weihenstephan_34_70_2,Sunv_A9,Suva_MCYC_623' \
        --skip Smik_IFO_1815_2='same strain and same taxon_id, keep one based on filtered sequence length' \
        --skip Spas_Weihenstephan_34_70_1='same strain and same taxon_id, keep one based on filtered sequence length' \
        --skip Skud_FM1057='short contigs' \
        --skip Skud_FM1062='short contigs' \
        --skip Skud_FM1066='short contigs' \
        --skip Skud_FM1069='short contigs' \
        --skip Skud_FM1071='short contigs' \
        --skip Skud_FM1072='short contigs' \
        --skip Skud_FM1073='short contigs' \
        --skip Skud_FM1074='short contigs' \
        --skip Skud_FM1075='short contigs' \
        --skip Skud_FM1076='short contigs' \
        --skip Skud_FM1076='short contigs' \
        --skip Skud_FM1077='short contigs' \
        --skip Skud_FM1078='short contigs' \
        --skip Skud_FM1079='short contigs' \
        --skip Skud_FM1094='short contigs' \
        --skip Skud_IFO10990='short contigs' \
        --skip Skud_IFO10991='short contigs' \
        --skip Skud_IFO1803='short contigs'
    ```

4. `pop_prep.pl`

    ```bash
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/saccharomyces_test.yml
    
    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh
    sh 04_plan_ALL.sh
    ```

5. `04_plan_ALL.sh`

    ```bash
    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

6. For other plans

    ```bash
    sh plan_XXX.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

7. Cleaning.
    
    * Remove useless files
    
        ```bash
        cd ~/data/alignment/saccharomyces
        
        # clean raw fasta
        find . -maxdepth 1 -type d -name "*_raw" | xargs rm -fr
        
        # clean maf-fasta
        find . -maxdepth 1 -type d -name "*_fasta" | xargs rm -fr
        
        # compress files
        find . -name "*.maf" | parallel gzip
        find . -name "*.fas" | parallel gzip
        ```

    * Restore everything to the beginning

        ```bash
        cd ~/data/alignment/saccharomyces

        find . -maxdepth 1 -type d -not -path "*WGS" | grep "\.\/" | xargs rm -fr
        rm *.xlsx *.csv *.sh *.bat *.nwk
        ```

## Candida WGS

1. Create pop/candida.tsv manually.

    http://www.ncbi.nlm.nih.gov/assembly?term=Candida

    Query a local ar_genbank DB.
    
    ```sql
    SELECT 
        CONCAT(LEFT(genus, 1),
                LEFT(TRIM(REPLACE(species, genus, '')),
                    3),
                REPLACE((REPLACE(organism_name, species, '')),
                    ' ',
                    '_')),
        SUBSTRING(wgs_master, 1, 4),
        organism_name,
        assembly_level
    FROM
        ar_genbank.ar
    WHERE
        genus = 'Candida'
            AND wgs_master LIKE '%000%'
    ORDER BY assembly_level , organism_name
    ```

2. Create working direcotry and download WGS sequences.

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

3. Download Candida dubliniensis CD36 and Candida orthopsilosis Co 90-125

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

3. 'gen_pop_conf.pl`

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

4. `pop_prep.pl`

    ```bash
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/candida_test.yml
    
    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh
    sh 04_plan_ALL.sh
    ```

5. `04_plan_ALL.sh`

    ```bash
    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 6_var_list.sh
    sh 7_multi_db_only.sh
    ```

6. For other plans

    ```bash
    sh plan_XXX.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

