## STEPS

1. Create `pop/trichoderma.tsv` manually. Names should only contain alphanumeric characters and underscores.
  Be careful with tabs and spaces, because .tsv stands for Tab-separated values, white spaces matters.

    Check NCBI pages

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Trichoderma&order=organism
    * http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5543
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid5543[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid5543[Organism:exp]

    And query a local ar_genbank DB. This is just a convenient but not accurate approach, especially for sub-species parts.

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
        genus = 'Trichoderma'
            AND wgs_master LIKE '%000%'
    ORDER BY assembly_level , organism_name
    ```

    When the two approach get very different number of strains, you can use the [downloaded list](http://www.ncbi.nlm.nih.gov/Traces/wgs/?&size=100&term=Trichoderma&state=live&retmode=text),
    create a quick .tsv file containing BioProject,Prefix,Organism,# of Contigs.

    ```bash
    curl 'http://www.ncbi.nlm.nih.gov/Traces/wgs/?&size=100&term=Trichoderma&retmode=text&size=all' \
        | perl -nl -a -F"\t" -e 'print qq{$F[2]\t$F[0]\t$F[4]\t$F[5]}' \
        > raw.tsv
    ```

    Then combined with the result from sql query. Run `wgs_prep.pl` to get `WGS/trichoderma.csv`,
    filter out duplicated items based on WGS prefix and run `wgs_prep.pl` again.

2. `wgs_prep.pl` will create a directory named `WGS` and three files containing meta information:

    ```bash
    mkdir -p ~/data/alignment/trichoderma
    cd ~/data/alignment/trichoderma

    perl ~/Scripts/withncbi/util/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/trichoderma.tsv \
        --fix \
        -o WGS \
        -a
    ```

    Working directory should be `~/data/alignment/trichoderma` throughout this procedure.

    1. `trichoderma.csv`

        Various information for WGS projects extracted from NCBI WGS record pages.

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

3. Download WGS files.

    ```bash
    # download with aria2
    aria2c -x 6 -s 3 -c -i WGS/trichoderma.url.txt

    # check downloaded .gz files
    find WGS -name "*.gz" | xargs gzip -t

    # rsync remote files
    # My connection to NCBI isn't stable, so download sequences in a linode VPS.
    # PLEASE don't hack it.
    # rsync --progress -av wangq@139.162.23.84:/home/wangq/data/alignment/trichoderma/ ~/data/alignment/trichoderma
    ```

4. Use `gen_pop_conf.pl` find matched files for each data entries in YAML and store extra options.

    * `per_seq` mean split fasta by names, target or good assembles should set it.
    * `skip` mean skip this strain.
    * `--opt rm_species=Fungi` specify the species or clade of this group for RepeatMasker.
    * `--opt min_contig=5000` to exclude short contigs.
    * `--opt per_seq_min_contig=20000` to exclude more short contigs for potential target.
    * The script will not overwrite existing .yml file by default, use `-y` to force it.

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
        --opt min_contig=5000 \
        --opt per_seq_min_contig=30000 \
        --plan 'name=four_way;t=Tatr_IMI_206040;qs=Tatr_XS215,Tree_QM6a,Tvir_Gv29_8' \
        --skip Tham_GD12='contigs are too short' \
        --per_seq Tatr_IMI_206040
    ```

5. Edit `pop/trichoderma_test.yml` on necessary.

    The following cmd refresh existing YAML file.

    ```bash
    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/Scripts/withncbi/pop/trichoderma_test.yml
    ```

6. `pop_prep.pl` will generate four or more bash scripts:

    ```bash
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/trichoderma_test.yml
    ```

    1. `01_file.sh`: unzip, filter and split
    2. `02_rm.sh`: RepeatMasker
    3. `03_strain_info.sh`: strain_info
    4. `plan_ALL.sh`: alignment plan for all genomes
    5. `plan_four_way.sh`: alignment plan for `four_way`, specified by `--plan` of `gen_pop_conf.pl`

7. Run generated scripts.

    Scripts starting with a "0" means they are doing preparing works.

    ```bash
    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh
    ```

8. `plan_ALL.sh` generates some bash files, execute them in order.

    ```bash
    sh plan_ALL.sh

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

10. When you are satisfied and don't see any wrong, rename `pop/trichoderma_test.yml` to
    `pop/trichoderma_data.yml` and commit it.

    ```bash
    mv ~/Scripts/withncbi/pop/trichoderma_test.yml ~/Scripts/withncbi/pop/trichoderma_data.yml
    ```

99. Cleaning. This is the ultimate final step. No more.

    * Remove useless files

    ```bash
    cd ~/data/alignment/trichoderma

    # clean raw fasta
    find . -maxdepth 1 -type d -name "*_raw" | xargs rm -fr

    # clean maf-fasta
    find . -maxdepth 1 -type d -name "*_fasta" | xargs rm -fr

    # clean raxml phy
    find . -maxdepth 2 -type f -name "*.phy" -or -name "*.phy.reduced" | xargs rm

    # compress files
    find . -type f -name "*.maf" | parallel gzip
    find . -type f -name "*.fas" | parallel gzip
    ```

    * Restore everything to the beginning

    ```bash
    cd ~/data/alignment/trichoderma

    find . -maxdepth 1 -type d -not -path "*WGS" | grep "\.\/" | xargs rm -fr
    rm *.xlsx *.csv *.sh *.bat *.nwk *.yml
    ```

## FAQ

* Why .tsv? All of your other programs use .csv.

    There are strains of which sub-species parts of names contain commas, can you believe it?

* I've 500 genomes of *E. coli*, the manually editing step 1 kills me.

    The whole `pop/` and much of `util/` scripts are for Eukaryotes. For small genomes
    of bacteria and archea, check `taxon/bac_prepare.pl`.
