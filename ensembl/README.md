# Ensembl related scripts

[TOC levels=1-3]: # " "

- [Ensembl related scripts](#ensembl-related-scripts)
  - [Download](#download)
    - [Ensembl mysql](#ensembl-mysql)
    - [EG mysql](#eg-mysql)
    - [Fasta](#fasta)
    - [Gff3](#gff3)
  - [Build local databases](#build-local-databases)
  - [Configurations](#configurations)
  - [`egaz prepseq`](#egaz-prepseq)


The current version of ensembl is 105 (Dec 2021) and the version of ensemblgenomes is 52.

## Download

Ensembl and Ensembl Genomes provides rsync service:

* `rsync://ftp.ensembl.org/ensembl/pub/release-105/`

* `rsync://ftp.ensemblgenomes.org/all/pub/release-52/plants/`

### Ensembl mysql

```bash
mkdir -p ~/data/ensembl105/mysql

for n in \
    homo_sapiens_core* \
    mus_musculus_core* \
    caenorhabditis_elegans_core* \
    drosophila_melanogaster_core* \
    saccharomyces_cerevisiae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        --exclude='*_cdna_*' \
        --exclude='*_vega_*' \
        --exclude='*_funcgen_*' \
        --exclude='*_otherfeatures_*' \
        --exclude='*_variation_*' \
        --exclude='*_rnaseq_*' \
        --exclude='ensembl_*' \
        --exclude='*_mart_*' \
        rsync://ftp.ensembl.org/ensembl/pub/release-105/mysql/${n} \
        ~/data/ensembl105/mysql
done

```

### EG mysql

```bash
mkdir -p ~/data/ensembl105/mysql

# Plants
for n in \
    arabidopsis_thaliana_core* arabidopsis_halleri_core* arabidopsis_lyrata_core* \
    brassica_napus_core* brassica_oleracea_core* brassica_rapa_core* \
    oryza_sativa_core* oryza_indica_core* \
    oryza_nivara_core* oryza_rufipogon_core* \
    solanum_lycopersicum_core* solanum_tuberosum_core* \
    vigna_angularis_core* vigna_radiata_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-52/plants/mysql/${n} \
        ~/data/ensembl105/mysql
done

# Metazoa
for n in \
    drosophila_sechellia_core* drosophila_simulans_core* \
    caenorhabditis_briggsae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-52/metazoa/mysql/${n} \
        ~/data/ensembl105/mysql
done

# Fungi
for n in \
    schizosaccharomyces_pombe_core* schizosaccharomyces_cryophilus_core* \
    schizosaccharomyces_japonicus_core* schizosaccharomyces_octosporus_core* \
    aspergillus_fumigatus_core* aspergillus_oryzae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-52/fungi/mysql/${n} \
        ~/data/ensembl105/mysql
done

# Protists
for n in \
    plasmodium_falciparum_core* dictyostelium_discoideum_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-52/protists/mysql/${n} \
        ~/data/ensembl105/mysql
done

```

### Fasta

```bash
mkdir -p ~/data/ensembl105/fasta

rsync -avP \
    --exclude='*ancestral_*' \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/release-105/fasta/ \
    ~/data/ensembl105/fasta

for n in \
    plants metazoa fungi protists \
    ; do
    echo "==> ${n}"
    rsync -avP \
        --exclude='caenorhabditis_elegans' \
        --exclude='drosophila_melanogaster' \
        --exclude='saccharomyces_cerevisiae' \
        --exclude='*ancestral_*' \
        --exclude='*_collection' \
        --exclude='*.dna.*' \
        --exclude='*.dna_rm.*' \
        --exclude='*.chromosome.*' \
        --exclude='*.nonchromosomal.*' \
        rsync://ftp.ensemblgenomes.org/all/pub/release-52/${n}/fasta/ \
        ~/data/ensembl105/fasta
done

```

### Gff3

```bash
mkdir -p ~/data/ensembl105/gff3

rsync -avP \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/release-105/gff3/ \
    ~/data/ensembl105/gff3

for n in \
    plants metazoa fungi protists \
    ; do
    echo "==> ${n}"
    rsync -avP \
        --exclude='caenorhabditis_elegans' \
        --exclude='drosophila_melanogaster' \
        --exclude='saccharomyces_cerevisiae' \
        --exclude='*_collection' \
        --exclude='*.dna.*' \
        --exclude='*.dna_rm.*' \
        --exclude='*.chromosome.*' \
        --exclude='*.nonchromosomal.*' \
        rsync://ftp.ensemblgenomes.org/all/pub/release-52/${n}/gff3/ \
        ~/data/ensembl105/gff3
done

```

## Build local databases

Use `build_ensembl.pl`.

```bash
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --checksum --ensembl ~/data/ensembl105/mysql/homo_sapiens_core_105_37

perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --initdb --db human_105 --ensembl ~/data/ensembl105/mysql/homo_sapiens_core_105_37

```

Or use `ensembl_batch.pl`, see [this](README.md#configurations).

## Configurations

Configurations stored in `ensembl_105.yml`.

```bash
cd ~/Scripts/withncbi/ensembl/
perl ensembl_batch.pl -i ensembl_105.yml

bash ensembl.build.sh
bash ensembl.fasta.sh
bash ensembl.anno.sh

cp ensembl.initrc.pm ~/Scripts/alignDB/

```

## `egaz prepseq`

```bash
find ~/data/alignment/Ensembl -maxdepth 1 -type d |
    sort |
    parallel --no-run-if-empty --linebuffer -k -j 2 '
        echo "==> {}"

        if [ -e {}/chr.2bit ]; then
            echo "    chr.2bit exists"
            exit
        fi

        if compgen -G "{}/*.fa" > /dev/null; then
            egaz prepseq {}
        else
            echo "    Can not find *.fa"
            exit
        fi
    '

```

