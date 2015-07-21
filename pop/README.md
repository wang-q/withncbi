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

5. Use 'gen_pop_conf.pl` find matched files for each @data entry in YAML and store extra options.

    * `per_seq` mean split fasta by names, target or good assembles should set it.
    * `skip` mean skip this strain.

    ```bash
    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/trichoderma/WGS/trichoderma.data.yml \
        -o ~/Scripts/withncbi/pop/trichoderma_test.yml \
        -d ~/data/alignment/trichoderma/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=trichoderma \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/trichoderma' \
        --opt rm_species=Fungi \
        --opt min_contig=10000 \
        --plan 'name=four_way;t=Tart_IMI_2206040;qs=Tatr_XS215,Tree_QM6a,Tvir_Gv29_8' \
        --skip Tham_GD12='contigs are too short' \
        --per_seq Tatr_IMI_2206040
    ```

6. Edit pop/trichoderma_test.yml on necessary.

    The following cmd refresh existing YAML file.

    ```bash
    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/Scripts/withncbi/pop/trichoderma_test.yml
    ```

7. `pop_prep.pl` will generate four or more bash scripts:

    ```bash
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/trichoderma_test.yml
    ```
    
    1. `01_file.sh`: unzip, filter and split
    2. `02_rm.sh`: RepeatMasker
    3. `03_strain_info.sh`: strain_info
    4. `04_plan_ALL.sh`: alignment plan for all genomes
    5. `plan_four_way.sh`: alignment plan for `four_way`

8. Run scripts.

    ```bash
    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh
    sh 04_plan_ALL.sh
    ```

8. `04_plan_ALL.sh` generated the following bash files, execute them.

    ```bash
    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

9. `plan_four_way.sh` will overwrite some bash files, execute the following:

    ```bash
    sh plan_four_way.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

10. When you are satisfied and don't see any wrong, rename pop/trichoderma_test.yml to
    pop/trichoderma_data.yml and git commit it.
    
    ```bash
    mv ~/Scripts/withncbi/pop/trichoderma_test.yml ~/Scripts/withncbi/pop/trichoderma_data.yml
    ```
