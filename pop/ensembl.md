# Downlaod Ensembl data

Current version of ensembl is 82, ensembl genomes' one is 29.

Ensembl provides rsync service.

Ensembl genomes only provide ftp.

## Ensembl

### mysql

```bash
mkdir -p ~/data/ensembl82/mysql
cd ~/data/ensembl82/mysql

rsync -avP \
    --exclude='*_cdna_82*' \
    --exclude='*_vega_82*' \
    --exclude='*_funcgen_82*' \
    --exclude='*_otherfeatures_82*' \
    --exclude='*_variation_82*' \
    --exclude='*_rnaseq_82*' \
    --exclude='ensembl_*' \
    --exclude='*_mart_82' \
    rsync://ftp.ensembl.org/ensembl/pub/release-82/mysql/ \
    .

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/grch37/release-82/mysql/homo_sapiens_core_82_37 .
```

### fasta

```bash
mkdir -p ~/data/ensembl82/fasta
cd ~/data/ensembl82/fasta

rsync -avP \
    --exclude='caenorhabditis_elegans' \
    --exclude='drosophila_melanogaster' \
    --exclude='saccharomyces_cerevisiae' \
    --exclude='homo_sapiens' \
    rsync://ftp.ensembl.org/ensembl/pub/release-82/fasta/ \
    .

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-75/fasta/homo_sapiens . # GRCh37
```

## Plants and Metazoa

### mysql

```bash
mkdir -p ~/data/ensembl82/mysql
cd ~/data/ensembl82/mysql

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/arabidopsis_thaliana_core_29_82_10 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/arabidopsis_lyrata_core_29_82_10 .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/oryza_sativa_core_29_82_7 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/oryza_indica_core_29_82_2 .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/brassica_oleracea_core_29_82_1 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/brassica_rapa_core_29_82_1 .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/solanum_lycopersicum_core_29_82_250 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/solanum_tuberosum_core_29_82_4 .

mv ftp.ensemblgenomes.org/pub/plants/release-29/mysql/* .
rm -fr ftp.ensemblgenomes.org
find . -name ".listing" | xargs rm

wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/drosophila_melanogaster_core_29_82_6 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/drosophila_sechellia_core_29_82_1 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/drosophila_simulans_core_29_82_1 .

wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_elegans_core_29_82_245 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_briggsae_core_29_82_230 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_remanei_core_29_82_233 .

mv ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/* .
rm -fr ftp.ensemblgenomes.org
find . -name ".listing" | xargs rm

# rsync -avP wangq@45.79.80.100:data/ensembl82/mysql/ ~/data/ensembl82/mysql
```

### fasta

```bash
mkdir -p ~/data/ensembl82/fasta
cd ~/data/ensembl82/fasta

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/arabidopsis_thaliana .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/arabidopsis_lyrata .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/oryza_sativa .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/oryza_indica .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/brassica_oleracea .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/brassica_rapa .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/solanum_lycopersicum .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/fasta/solanum_tuberosum .

mv ftp.ensemblgenomes.org/pub/plants/release-29/fasta/* .
rm -fr ftp.ensemblgenomes.org
find . -name ".listing" | xargs rm

wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/fasta/drosophila_melanogaster .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/fasta/drosophila_sechellia .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/fasta/drosophila_simulans .

wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/fasta/caenorhabditis_elegans .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/fasta/caenorhabditis_briggsae .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/fasta/caenorhabditis_remanei .

mv ftp.ensemblgenomes.org/pub/metazoa/release-29/fasta/* .
rm -fr ftp.ensemblgenomes.org
find . -name ".listing" | xargs rm

# rsync -avP wangq@45.79.80.100:data/ensembl82/fasta/ ~/data/ensembl82/fasta
```
