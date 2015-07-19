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

5. Add some fields to WGS/trichoderma.data.yml, get pop/trichoderma_data.yml.

    * `per_seq` mean split fasta by names, target or good assembles should set it.
    * `skip` mean skip this strains.
    * The rest are running parameters.

    ```yaml
    ---
    data:
        - taxon: 452589
          name: 'Tart_IMI_2206040'
          per_seq: 1
        - taxon: 1247866
          name: 'Tham_GD12'
          skip: 'contigs are too short'
    group_name: 'trichoderma'
    base_dir: '~/data/alignment'
    data_dir: '~/data/alignment/trichoderma'
    fasta_dir: '~/data/alignment/trichoderma/WGS'
    pl_dir: '~/Scripts'
    parallel: 4
    ```

6. Edit trichoderma.pl

    * Section `01_file`: length thresholds
    * Section `02_rm`: -species Fungi
    * Section `03_prepare`: multi genome alignment plan

7. `trichoderma.pl` will match entries in yaml file and downloaded files, then generate three bash scripts:

    `perl ~/Scripts/withncbi/pop/trichoderma.pl`
    
    1. `01_file.sh`: unzip, filter and split
    2. `02_rm.sh`: RepeatMasker
    3. `03_prepare.sh`: aligning plans
        * copy & paste lines of `03_prepare.sh` section by section into terminal.

8. For each aligning plans (multi_name), execute the following bash file.

    ```bash
    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 6_var_list.sh
    sh 7_multi_db_only.sh
    sh 9_pack_it_up.sh
    ```
