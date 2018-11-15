# Ensembl related scripts.

Current version of ensembl is 94, the one of ensembl genomes is 41.

[TOC levels=1-3]: # " "
- [Ensembl related scripts.](#ensembl-related-scripts)
- [Downlaod Ensembl data](#downlaod-ensembl-data)
    - [Ensembl mysql](#ensembl-mysql)
    - [Ensembl fasta](#ensembl-fasta)
    - [Ensembl gff3](#ensembl-gff3)
- [Ensembl genomes](#ensembl-genomes)
    - [EG mysql](#eg-mysql)
    - [EG fasta](#eg-fasta)
    - [EG gff3](#eg-gff3)
- [Build local databases](#build-local-databases)
- [Configurations](#configurations)


# Downlaod Ensembl data

Ensembl and Ensembl Genomes provides rsync service:

* `rsync://ftp.ensembl.org/ensembl/pub/release-94/`

* `rsync://ftp.ensemblgenomes.org/all/pub/plants/release-41/`

## Ensembl mysql

```bash
mkdir -p ~/data/ensembl94/mysql

for n in \
    pan_troglodytes_core* gorilla_gorilla_core* \
    pongo_abelii_core* macaca_mulatta_core* \
    mus_musculus_core* rattus_norvegicus_core* \
    caenorhabditis_elegans_core* drosophila_melanogaster_core* \
    saccharomyces_cerevisiae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        --exclude='*_cdna_94*' \
        --exclude='*_vega_94*' \
        --exclude='*_funcgen_94*' \
        --exclude='*_otherfeatures_94*' \
        --exclude='*_variation_94*' \
        --exclude='*_rnaseq_94*' \
        --exclude='ensembl_*' \
        --exclude='*_mart_94' \
        --exclude='homo_sapiens_core_94_38' \
        rsync://ftp.ensembl.org/ensembl/pub/release-94/mysql/${n} \
        ~/data/ensembl94/mysql
done

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/grch37/release-94/mysql/homo_sapiens_core_94_37 ~/data/ensembl94/mysql

```

## Ensembl fasta

```bash
mkdir -p ~/data/ensembl94/fasta

rsync -avP \
    --exclude='homo_sapiens' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/release-94/fasta/ \
    ~/data/ensembl94/fasta

# GRCh37
rsync -avP \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/grch37/release-94/fasta/homo_sapiens \
    ~/data/ensembl94/fasta

```

## Ensembl gff3

```bash
mkdir -p ~/data/ensembl94/gff3

rsync -avP \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/release-94/gff3/ \
    ~/data/ensembl94/gff3

```

# Ensembl genomes

##  EG mysql

```bash
mkdir -p ~/data/ensembl94/mysql

# Plants
for n in \
    arabidopsis_thaliana_core* arabidopsis_lyrata_core* \
    oryza_sativa_core* oryza_indica_core* \
    brassica_oleracea_core* brassica_rapa_core* \
    solanum_lycopersicum_core* solanum_tuberosum_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-41/plants/mysql/${n} \
        ~/data/ensembl94/mysql
done

# Metazoa
for n in \
    drosophila_sechellia_core* drosophila_simulans_core* \
    caenorhabditis_briggsae_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-41/metazoa/mysql/${n} \
        ~/data/ensembl94/mysql
done

# Fungi
for n in \
    schizosaccharomyces_pombe_core* aspergillus_fumigatus_core* \
    ensembl_compara_fungi_41_94 \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-41/fungi/mysql/${n} \
        ~/data/ensembl94/mysql
done

# Protists
for n in \
    plasmodium_falciparum_core* dictyostelium_discoideum_core* \
    ; do
    echo "==> ${n}"
    rsync -avP \
        rsync://ftp.ensemblgenomes.org/all/pub/release-41/protists/mysql/${n} \
        ~/data/ensembl94/mysql
done

# compara
rsync -avP \
    rsync://ftp.ensemblgenomes.org/all/pub/release-41/fungi/mysql/ensembl_compara_fungi_41_94 \
    ~/data/ensembl94/mysql

```

## EG fasta

```bash
mkdir -p ~/data/ensembl94/fasta

for n in \
    plants metazoa fungi protists \
    ; do
    rsync -avP \
        --exclude='caenorhabditis_elegans' \
        --exclude='drosophila_melanogaster' \
        --exclude='saccharomyces_cerevisiae' \
        --exclude='*_collection' \
        --exclude='*.dna.*' \
        --exclude='*.dna_rm.*' \
        --exclude='*.chromosome.*' \
        --exclude='*.nonchromosomal.*' \
        rsync://ftp.ensemblgenomes.org/all/pub/release-41/${n}/fasta/ \
        ~/data/ensembl94/fasta
done

```

## EG gff3

```bash
mkdir -p ~/data/ensembl94/gff3

for n in \
    plants metazoa fungi protists \
    ; do
    rsync -avP \
        --exclude='caenorhabditis_elegans' \
        --exclude='drosophila_melanogaster' \
        --exclude='saccharomyces_cerevisiae' \
        --exclude='*_collection' \
        --exclude='*.dna.*' \
        --exclude='*.dna_rm.*' \
        --exclude='*.chromosome.*' \
        --exclude='*.nonchromosomal.*' \
        rsync://ftp.ensemblgenomes.org/all/pub/release-41/${n}/gff3/ \
        ~/data/ensembl94/gff3
done

```

# Build local databases

Use `build_ensembl.pl`.

```bash
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --checksum --ensembl ~/data/ensembl94/mysql/homo_sapiens_core_94_37

perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --initdb --db human_94 --ensembl ~/data/ensembl94/mysql/homo_sapiens_core_94_37
```

Or use `ensembl_batch.pl`, see [this](README.md#configurations).

# Configurations

Configurations stored in `ensembl_94.yml`.

```bash
cd ~/Scripts/withncbi/ensembl/
perl ensembl_batch.pl -i ensembl_94.yml

bash ensembl.build.sh
bash ensembl.fasta.sh
bash ensembl.anno.sh

cp ensembl.initrc.pm ~/Scripts/alignDB/
```

