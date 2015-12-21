# Build alignments of several strains of *Saccharomyces cerevisiae*

As example for egaz and alignDB. Extract from Scer_wgs of `OP-*.md`.

## Download

1. Create `scer_example.tsv` manually.

    ```bash
    mkdir -p ~/data/alignment/example/scer          # operation directory
    mkdir -p ~/data/alignment/example/GENOMES  # sequence directory

    cd ~/data/alignment/example/GENOMES

    echo -e '#name\tprefix\torganism\tcontigs' > example_wgs.tsv
    echo -e "Spar\tAABY\tSaccharomyces paradoxus NRRL Y-17217\t832" >> example_wgs.tsv
    ```

2. Create working directory and download WGS sequences.

    ```bash
    cd ~/data/alignment/example/GENOMES

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f example_wgs.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/example_wgs.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download strains of *Saccharomyces cerevisiae* at good assembly status.

    ```bash
    mkdir -p ~/data/alignment/example/GENOMES/DOWNLOAD
    cd ~/data/alignment/example/GENOMES/DOWNLOAD

    # Download S288c, EC1118 and RM11_1a separately
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt \
        --nuclear -name S288c \
        > S288c.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000218975.1.assembly.txt \
        --nuclear --genbank --scaffold -name EC1118 \
        > EC1118.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000149365.1.assembly.txt \
        --genbank --scaffold -name RM11_1a \
        > RM11_1a.seq.csv

    echo "#strain_name,accession,strain_taxon_id,seq_name" > example.seq.csv
    cat S288c.seq.csv EC1118.seq.csv RM11_1a.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> example.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f example.seq.csv

    ```

## Align

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/example/scer
    cd ~/data/alignment/example/scer

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/example/GENOMES/WGS/example_wgs.data.yml \
        -o scer_test.yml \
        -d ~/data/alignment/example/GENOMES/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=scer \
        --opt base_dir='~/data/alignment/example' \
        --opt data_dir="~/data/alignment/example/scer" \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/example/GENOMES/DOWNLOAD \
        --download "name=S288c;taxon=559292" \
        --download "name=RM11_1a;taxon=285006" \
        --download "name=EC1118;taxon=643680" \
        --plan 'name=Scer_n3_Spar;t=S288c;qs=EC1118,RM11_1a,Spar;o=Spar' \
        --plan 'name=Scer_n3_pop;t=S288c;qs=EC1118,RM11_1a' \
        -y
    ```

2. Rest routing things.

    ```bash
    cd ~/data/alignment/example/scer

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i scer_test.yml

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
    sh plan_Scer_n3_Spar.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Scer_n3_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```
