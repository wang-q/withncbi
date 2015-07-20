## STEPS

1. Create pop/trichoderma.tsv manually, be careful with tabs and spaces.

    http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=trichoderma

2. Working directory should be `~/data/alignment/trichoderma` throughout this procedure.

    ```bash
    mkdir -p ~/data/alignment/trichoderma
    cd ~/data/alignment/trichoderma
    ```

3. `wgs_prep.pl` will create a directory named `WGS` and three files:

    ```bash
    perl ~/Scripts/withncbi/util/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/trichoderma.tsv \
        --fix \
        -o WGS \
        -a
    ```

    1. `trichoderma.csv`
    
        Variaous infomation for WGS projects extracted from NCBI WGS record pages.
    
    2. `trichoderma.data.yml`
    
        ```yaml
        ---
        data:
            - taxon: 452589
              name: 'Tart_IMI_2206040'
              sciname: 'Trichoderma atroviride IMI 206040'
              prefix: 'ABDG02'
              coverage: '8.26x Sanger'
        ```
        
    3. `trichoderma.url.txt`
    
        Download urls for WGS files.

4. Download WGS files.

    ```bash
    # download with aria2
    aria2c -x 6 -s 3 -c -i WGS/trichoderma.url.txt
    
    # check downloaded .gz files
    find WGS -name "*.gz" | xargs gzip -t 
    ```

5. Use 'match_data.pl` find matched files for each @data entry in YAML and store extra options.

    * `per_seq` mean split fasta by names, target or good assembles should set it.
    * `skip` mean skip this strain.

    ```bash
    perl ~/Scripts/withncbi/pop/match_data.pl \
        -i ~/data/alignment/trichoderma/WGS/trichoderma.data.yml \
        -o ~/Scripts/withncbi/pop/trichoderma_test.yml \
        -d ~/data/alignment/trichoderma/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=trichoderma \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/trichoderma' \
        --opt rm_species=Fungi \
        --skip Tham_GD12='contigs are too short' \
        --per_seq Tart_IMI_2206040
    ```

7. Add multiply alignment plans to pop/trichoderma_test.yml.

    ```yaml
    ---
    plans:
        - name: 'four_way'
          t: 'Tart_IMI_2206040'
          qs:
            - 'Tvir_Gv29_8'

    ```


8. `pop_prep.pl` will generate three or more bash scripts:

    `perl ~/Scripts/withncbi/pop/pop_prep.pl -i trichoderma_test.yml`
    
    1. `01_file.sh`: unzip, filter and split
    2. `02_rm.sh`: RepeatMasker
    3. `03_strain_info.sh`: strain_info and alignment plan of all genomes
    4. `plan_four_way.sh`: alignment plan of `four_way`

9. For each aligning plans (multi_name), execute the following bash file.

    ```bash
    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 6_var_list.sh
    sh 7_multi_db_only.sh
    sh 9_pack_it_up.sh
    ```
