# Operating steps for each group.

## Candida

1. Create pop/candida.tsv manually, be careful with tabs and spaces.

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
            AND assembly_level != 'Contig'
            AND wgs_master LIKE '%000%'
    ORDER BY organism_name
    ```
2. Create working direcotry and download sequences.

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
3. 'match_data.pl`

    ```bash
    perl ~/Scripts/withncbi/pop/match_data.pl \
        -i ~/data/alignment/candida/WGS/candida.data.yml \
        -o ~/Scripts/withncbi/pop/candida_test.yml \
        -d ~/data/alignment/candida/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=candida \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/candida' \
        --opt rm_species=Fungi \
        --per_seq Calb_WO_1
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
