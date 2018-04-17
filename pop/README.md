# Build alignments on an whole Eukaryotes genus

Or order, family, species.

Genus Trichoderma as example.

[TOC levels=1-3]: # " "
- [Build alignments on an whole Eukaryotes genus](#build-alignments-on-an-whole-eukaryotes-genus)
- [Section 1: select strains and download sequences.](#section-1-select-strains-and-download-sequences)
    - [`pop/trichoderma.tsv`](#poptrichodermatsv)
    - [`wgs_prep.pl`](#wgs_preppl)
    - [Download WGS files.](#download-wgs-files)
    - [Download ASSEMBLY files](#download-assembly-files)
- [Section 2: create configuration file and generate alignments.](#section-2-create-configuration-file-and-generate-alignments)
- [Section 3: cleaning.](#section-3-cleaning)
- [FAQ](#faq)


# Section 1: select strains and download sequences.

## `pop/trichoderma.tsv`

Create `pop/trichoderma.tsv` manually. Names should only contain alphanumeric characters and
underscores. Be careful with tabs and spaces, because .tsv stands for Tab-separated values, white
spaces matters.

Check NCBI pages

* https://www.ncbi.nlm.nih.gov/Traces/wgs/?view=wgs&search=Trichoderma
* http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5543
* http://www.ncbi.nlm.nih.gov/genome/?term=txid5543[Organism:exp]
* http://www.ncbi.nlm.nih.gov/assembly?term=txid5543[Organism:exp]

And query a local `ar_genbank` DB. This is just a convenient but not accurate approach, especially
for sub-species parts.

```mysql
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

For genus contains many species, you should be careful that "Gspe" (*G*enus *spe*cies) style
abbreviation may mix up two or more species.

```mysql
SELECT DISTINCT
    species
FROM
    ar_genbank.ar
WHERE
    genus = 'Trichoderma'
```

When the two approaches get very different number of strains, you run the following steps. Check
intermediate results on necessary.

Working directory should be `~/data/alignment/trichoderma` in this section.

```bash
export RANK_LEVEL=genus
export RANK_ID=5543
export RANK_NAME=trichoderma

mkdir -p ~/data/alignment/${RANK_NAME}            # Working directory

cd ~/data/alignment/${RANK_NAME}
```

You can copy & paste the following block of codes as a whole unit.

```bash
# stage1
# Results from sql query.
mysql -ualignDB -palignDB gr_euk -e "
    SELECT 
        SUBSTRING(wgs,1,6) as prefix0,
        SUBSTRING(wgs,1,4) as prefix,
        organism_name,
        status 
    FROM gr 
    WHERE wgs != '' AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        CONCAT(SUBSTRING(wgs_master,1,5), RIGHT(wgs_master,1)) as prefix0,
        SUBSTRING(wgs_master,1,4) as prefix,
        organism_name,
        assembly_level 
    FROM ar 
    WHERE wgs_master != '' AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    >> raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        CONCAT(SUBSTRING(wgs_master,1,5), RIGHT(wgs_master,1)) as prefix0,
        SUBSTRING(wgs_master,1,4) as prefix,
        organism_name,
        assembly_level 
    FROM ar 
    WHERE wgs_master != '' AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    >> raw.tsv

# stage2
# Click the 'Download' button in the middle of NCBI WGS page.
# NCBI changed its WGS page, make curl defunct
rm -f ~/Downloads/wgs_selector*.csv
chrome "https://www.ncbi.nlm.nih.gov/Traces/wgs/?view=wgs&search=${RANK_NAME}"

# Quit chrome Cmd-Q

# There're no chromosome level assemblies in WGS
cat ~/Downloads/wgs_selector.csv |
    perl -nl -a -F"," -e '
        my $p = substr($F[0],0,4);
        my $status = $F[13] > 0 ? q{Scaffold} : q{Contig};
        print qq{$F[0]\t$p\t$F[4]\t$status}
    ' \
    >> raw.tsv

cat raw.tsv |
    perl -nl -a -F"\t" -e '
        BEGIN{my %seen}; 
        /^prefix/i and next;
        scalar @F == 4 or next;
        $seen{$F[1]}++;
        $seen{$F[1]} > 1 and next;
        print join(qq{\t}, $F[0], $F[0], $F[2], $F[3]);
    ' \
    > raw2.tsv

# stage3
# Run `wgs_prep.pl` to get a crude `raw2.csv`
perl ~/Scripts/withncbi/taxon/wgs_prep.pl -f raw2.tsv --csvonly

echo -e '#name\tprefix\torganism\tcontigs' > raw3.tsv
cat raw2.csv |
    perl -nl -a -F"," -e '
        /^prefix/i and next;
        s/"//g for @F;
        @O = split(/ /, $F[3]);
        $F[4] =~ s/\s+$//g;
        $F[4] =~ s/\W+/_/g;
        $name = substr($O[0],0,1) . substr($O[1],0,3);
        $name .= q{_} . $F[4] if $F[4];
        print qq{$name\t$F[0]\t$F[3]\t$F[9]}
    ' |
    sort -t$'\t' -k4 -n |
    uniq \
    >> raw3.tsv

mv raw3.tsv ${RANK_NAME}.wgs.tsv

# find potential duplicated strains or assemblies
cat ${RANK_NAME}.wgs.tsv |
    perl -nl -a -F"\t" -e 'print $F[0]' |
    uniq -c

# Edit .tsv, remove duplicated strains, check strain names and comment out poor assemblies.
# vim ${GENUS}.wgs.tsv

```

Put the .tsv file to `~/Scripts/withncbi/pop/` and run `wgs_prep.pl` again. When everything is fine,
commit the .tsv file.

For detailed WGS info, click Prefix column lead to WGS project, where we could download gzipped
fasta and project description manually.

## `wgs_prep.pl`

`wgs_prep.pl` will create a directory named `WGS` and three files containing meta information:

```bash
cd ~/data/alignment/trichoderma

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/withncbi/pop/trichoderma.tsv \
    --fix \
    -o WGS \
    -a

```

1. `trichoderma.csv`

    Various information for WGS projects extracted from NCBI WGS record pages.

2. `trichoderma.url.txt`

    Download urls for WGS files.

3. `trichoderma.data.yml`

    ```yaml
    ---
    data:
      - taxon: 452589
        name: 'Tart_IMI_2206040'
        sciname: 'Trichoderma atroviride IMI 206040'
        prefix: 'ABDG02'
        coverage: '8.26x Sanger'
    ```

## Download WGS files.

```bash
# download with aria2
cd ~/data/alignment/trichoderma
aria2c -UWget -x 6 -s 3 -c -i WGS/trichoderma.url.txt

# check downloaded .gz files
find WGS -name "*.gz" | parallel -j 4 gzip -t

# rsync remote files
# My connection to NCBI isn't stable, so download sequences in a linode VPS.
# PLEASE don't hack it.
# rsync -avP wangq@173.230.144.105:data/alignment/trichoderma/ ~/data/alignment/trichoderma
```

## Download ASSEMBLY files


```bash

echo -e '#name\tftp_path\torganism\tassembly_level' > ${RANK_NAME}.assembly.tsv

# comment out unneeded conditions
mysql -ualignDB -palignDB ar_genbank -e "
    SELECT 
        organism_name, species, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
#        AND wgs_master = ''
#        AND assembly_level = 'Chromosome'
        AND organism_name != species
        AND ${RANK_LEVEL}_id = ${RANK_ID}
    " |
    perl -nl -a -F"\t" -e '
        /^organism_name/i and next;
        $n = $F[0];
        $rx = quotemeta $F[1];
        $n =~ s/$rx\s*//;
        $n =~ s/\s+$//;
        $n =~ s/\W+/_/g;
        @O = split(/ /, $F[1]);
        $name = substr($O[0],0,1) . substr($O[1],0,3);
        $name .= q{_} . $n if $n;
        printf qq{%s\t%s\t%s\t%s\n}, $name, $F[2], $F[1], $F[3];
        ' \
    >> ${RANK_NAME}.assembly.tsv

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ${RANK_NAME}.assembly.tsv \
    -o ASSEMBLY

```

```bash
cd ~/data/alignment/trichoderma

bash ASSEMBLY/trichoderma.assembly.rsync.sh

# rsync -avP wangq@173.230.144.105:data/alignment/trichoderma/ ~/data/alignment/trichoderma

bash ASSEMBLY/trichoderma.assembly.collect.sh

```

```bash
# Cleaning
rm raw*.*sv
unset RANK_LEVEL
unset RANK_ID
unset RANK_NAME

```

# Section 2: create configuration file and generate alignments.

Working directory should be `~/data/alignment/Fungi/trichoderma` in this section.

1. Use `gen_pop_conf.pl` find matched files for each data entries in YAML and store extra options.

    * `per_seq` mean split fasta by names, target or good assembles should set it.
    * `skip` mean skip this strain.
    * `--opt rm_species=Fungi` specify the species or clade of this group for RepeatMasker.
    * `--opt min_contig=5000` to exclude short contigs.
    * `--opt per_seq_min_contig=20000` to exclude more short contigs for potential target.
    * The script will not overwrite existing .yml file by default, use `-y` to force it.

    ```bash
    mkdir -p ~/data/alignment/Fungi/trichoderma
    cd ~/data/alignment/Fungi/trichoderma

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/GENOMES/trichoderma/WGS/trichoderma.data.yml \
        -o ~/Scripts/withncbi/pop/trichoderma_test.yml \
        -d ~/data/alignment/Fungi/GENOMES/trichoderma/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=trichoderma \
        --opt base_dir='~/data/alignment/Fungi/' \
        --opt data_dir='~/data/alignment/Fungi/trichoderma' \
        --opt rm_species=Fungi \
        --opt min_contig=5000 \
        --opt per_seq_min_contig=30000 \
        --plan 'name=four_way;t=Tatr_IMI_206040;qs=Tatr_XS2015,Tree_QM6a,Tvir_Gv29_8' \
        --skip Tham_GD12='contigs are too short' \
        --per_seq Tatr_IMI_206040
    ```

2. Edit `pop/trichoderma_test.yml` on necessary.

    The following cmd refresh existing YAML file.

    ```bash
    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/Scripts/withncbi/pop/trichoderma_test.yml
    ```

3. `pop_prep.pl` will generate four or more bash scripts:

    ```bash
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/trichoderma_test.yml
    ```

    1. `01_file.sh`: unzip, filter and split
    2. `02_rm.sh`: RepeatMasker
    3. `03_strain_info.sh`: strain_info
    4. `plan_ALL.sh`: alignment plan for all genomes
    5. `plan_four_way.sh`: alignment plan for `four_way`, specified by `--plan` of `gen_pop_conf.pl`

4. Run generated scripts.

    Scripts starting with a "0" means they are doing preparing works.

    ```bash
    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh
    ```

5. `plan_ALL.sh` generates some bash files, execute them in order.

    ```bash
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

6. `plan_four_way.sh` will overwrite some bash files, execute the following:

    ```bash
    sh plan_four_way.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

7. When you are satisfied and don't see any wrong, rename `pop/trichoderma_test.yml` to
   `pop/trichoderma_data.yml` and commit it.

    ```bash
    mv ~/Scripts/withncbi/pop/trichoderma_test.yml ~/Scripts/withncbi/pop/trichoderma_data.yml
    ```

# Section 3: cleaning.

This is the ultimate final step. No more. Actually you may choose not to do this. It's depended on
your disk capacity.

```bash
cd ~/data/alignment/Fungi/trichoderma

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

# FAQ

* Why .tsv? All of your other programs use .csv.

    There are strains of which sub-species parts contain commas, can you believe it?

* I've 500 genomes of *E. coli*, the manually editing step 1 kills me.

    The whole `pop/` and much of `util/` scripts are for Eukaryotes. For small genomes of bacteria,
    archaea and organelles, check `taxon/bac_prepare.pl`.

* Your command lines executed and the results are wired.

    Be sure you aren't in Windows. Be sure you are familiar to bash command lines.

    Or send what you want to me and let me do the job.

* I have a very good assembly on chromosome level, but I can't find it in WGS.

    Best genomes on the world went to NCBI RefSeq. Use tools in `util/` to download them. Examples
    can be found in `pop/OPs-download.md`.
