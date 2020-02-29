# Ensembl related scripts

Current version of ensembl is 98 (September 2019), the one of ensembl genomes is 45.

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


## Download

Ensembl and Ensembl Genomes provides rsync service:

* `rsync://ftp.ensembl.org/ensembl/pub/release-98/`

* `rsync://ftp.ensemblgenomes.org/all/pub/release-45/plants/`

### Ensembl mysql

```bash
mkdir -p ~/data/ensembl98/mysql

for n in \
    homo_sapiens_core* \
    pan_troglodytes_core* gorilla_gorilla_core* \
    pongo_abelii_core* macaca_mulatta_core* \
    mus_musculus_core* rattus_norvegicus_core* \
    caenorhabditis_elegans_core* drosophila_melanogaster_core* \
    saccharomyces_cerevisiae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        --exclude='*_cdna_98*' \
        --exclude='*_vega_98*' \
        --exclude='*_funcgen_98*' \
        --exclude='*_otherfeatures_98*' \
        --exclude='*_variation_98*' \
        --exclude='*_rnaseq_98*' \
        --exclude='ensembl_*' \
        --exclude='*_mart_98' \
        rsync://ftp.ensembl.org/ensembl/pub/release-98/mysql/${n} \
        ~/data/ensembl98/mysql
done

```

### EG mysql

```bash
mkdir -p ~/data/ensembl98/mysql

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
        rsync://ftp.ensemblgenomes.org/all/pub/release-45/plants/mysql/${n} \
        ~/data/ensembl98/mysql
done

# Metazoa
for n in \
    drosophila_sechellia_core* drosophila_simulans_core* \
    caenorhabditis_briggsae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-45/metazoa/mysql/${n} \
        ~/data/ensembl98/mysql
done

# Fungi
for n in \
    schizosaccharomyces_pombe_core* schizosaccharomyces_cryophilus_core* \
    schizosaccharomyces_japonicus_core* schizosaccharomyces_octosporus_core* \
    aspergillus_fumigatus_core* aspergillus_oryzae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-45/fungi/mysql/${n} \
        ~/data/ensembl98/mysql
done

# Protists
for n in \
    plasmodium_falciparum_core* dictyostelium_discoideum_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-45/protists/mysql/${n} \
        ~/data/ensembl98/mysql
done

```

### Fasta

```bash
mkdir -p ~/data/ensembl98/fasta

rsync -avP \
    --exclude='*ancestral_*' \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/release-98/fasta/ \
    ~/data/ensembl98/fasta

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
        rsync://ftp.ensemblgenomes.org/all/pub/release-45/${n}/fasta/ \
        ~/data/ensembl98/fasta
done

```

### Gff3

```bash
mkdir -p ~/data/ensembl98/gff3

rsync -avP \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/release-98/gff3/ \
    ~/data/ensembl98/gff3

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
        rsync://ftp.ensemblgenomes.org/all/pub/release-45/${n}/gff3/ \
        ~/data/ensembl98/gff3
done

```

## Build local databases

Use `build_ensembl.pl`.

```bash
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --checksum --ensembl ~/data/ensembl98/mysql/homo_sapiens_core_98_37

perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --initdb --db human_98 --ensembl ~/data/ensembl98/mysql/homo_sapiens_core_98_37

```

Or use `ensembl_batch.pl`, see [this](README.md#configurations).

## Configurations

Configurations stored in `ensembl_98.yml`.

```bash
cd ~/Scripts/withncbi/ensembl/
perl ensembl_batch.pl -i ensembl_98.yml

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

