# Ensembl related scripts.

Current version of ensembl is 82, ensembl genomes' one is 29.

Configurations stored in `ensembl_82.yml`.

```bash
cd ~/Scripts/withncbi/ensembl/
perl ensembl_batch.pl -i ensembl_82.yml

bash ensembl.build.sh
bash ensembl.fasta.sh
bash ensembl.anno.sh

cp ensembl.initrc.pm ~/Scripts/alignDB/
```

## Downlaod Ensembl data

Ensembl provides rsync service: `rsync://ftp.ensembl.org/ensembl/pub/release-82/`.

Ensembl genomes only provide ftp.

### Ensembl mysql

```bash
mkdir -p ~/data/ensembl82/mysql

rsync -avP \
    --exclude='*_cdna_82*' \
    --exclude='*_vega_82*' \
    --exclude='*_funcgen_82*' \
    --exclude='*_otherfeatures_82*' \
    --exclude='*_variation_82*' \
    --exclude='*_rnaseq_82*' \
    --exclude='ensembl_*' \
    --exclude='*_mart_82' \
    --exclude='homo_sapiens_core_82_38' \
    rsync://ftp.ensembl.org/ensembl/pub/release-82/mysql/ \
    ~/data/ensembl82/mysql

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/grch37/release-82/mysql/homo_sapiens_core_82_37 ~/data/ensembl82/mysql

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-82/mysql/ensembl_compara_82 ~/data/ensembl82/mysql
```

### Ensembl fasta

```bash
mkdir -p ~/data/ensembl82/fasta

rsync -avP \
    --exclude='caenorhabditis_elegans' \
    --exclude='drosophila_melanogaster' \
    --exclude='saccharomyces_cerevisiae' \
    --exclude='homo_sapiens' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensembl.org/ensembl/pub/release-82/fasta/ \
    ~/data/ensembl82/fasta

rsync -avP rsync://ftp.ensembl.org/ensembl/pub/release-75/fasta/homo_sapiens ~/data/ensembl82/fasta # GRCh37
```

### Ensembl genomes mysql

```bash
mkdir -p ~/data/ensembl82/mysql
cd ~/data/ensembl82/mysql

# Plants
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/arabidopsis_thaliana_core_29_82_10 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/arabidopsis_lyrata_core_29_82_10 .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/oryza_sativa_core_29_82_7 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/oryza_indica_core_29_82_2 .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/brassica_oleracea_core_29_82_1 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/brassica_rapa_core_29_82_1 .

wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/solanum_lycopersicum_core_29_82_250 .
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/mysql/solanum_tuberosum_core_29_82_4 .

# Metazoa
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/drosophila_melanogaster_core_29_82_6 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/drosophila_sechellia_core_29_82_1 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/drosophila_simulans_core_29_82_1 .

wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_elegans_core_29_82_245 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_briggsae_core_29_82_230 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_remanei_core_29_82_233 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_brenneri_core_29_82_233 .
wget -m ftp://ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/caenorhabditis_japonica_core_29_82_233 .

# Fungi
wget -m ftp://ftp.ensemblgenomes.org/pub/fungi/release-29/mysql/saccharomyces_cerevisiae_core_29_82_4 .
wget -m ftp://ftp.ensemblgenomes.org/pub/fungi/release-29/mysql/schizosaccharomyces_pombe_core_29_82_2 .
wget -m ftp://ftp.ensemblgenomes.org/pub/fungi/release-29/mysql/aspergillus_fumigatus_core_29_82_2 .

# Protists
wget -m ftp://ftp.ensemblgenomes.org/pub/protists/release-29/mysql/plasmodium_falciparum_core_29_82_3 .
wget -m ftp://ftp.ensemblgenomes.org/pub/protists/release-29/mysql/dictyostelium_discoideum_core_29_82_1 .

# compara
wget -m ftp://ftp.ensemblgenomes.org/pub/release-29/fungi/mysql/ensembl_compara_fungi_29_82 .

# clean
mv ftp.ensemblgenomes.org/pub/metazoa/release-29/mysql/* .
mv ftp.ensemblgenomes.org/pub/plants/release-29/mysql/* .
mv ftp.ensemblgenomes.org/pub/fungi/release-29/mysql/* .
mv ftp.ensemblgenomes.org/pub/protists/release-29/mysql/* .

find . -name ".listing" | xargs rm

# rsync -avP wangq@45.79.80.100:data/ensembl82/mysql/ ~/data/ensembl82/mysql
```

### Ensembl genomes fasta

```bash
mkdir -p ~/data/ensembl82/fasta
cd ~/data/ensembl82/fasta

# Plants
rsync -avP \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensemblgenomes.org/all/pub/plants/release-29/fasta/ \
    ~/data/ensembl82/fasta

# Protists
rsync -avP \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensemblgenomes.org/all/pub/protists/release-29/fasta/ \
    ~/data/ensembl82/fasta

# Fungi
rsync -avP \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensemblgenomes.org/all/pub/fungi/release-29/fasta/ \
    ~/data/ensembl82/fasta

# Metazoa
rsync -avP \
    --exclude='*_collection' \
    --exclude='*.dna.*' \
    --exclude='*.dna_rm.*' \
    --exclude='*.chromosome.*' \
    --exclude='*.nonchromosomal.*' \
    rsync://ftp.ensemblgenomes.org/all/pub/metazoa/release-29/fasta/ \
    ~/data/ensembl82/fasta

# rsync -avP wangq@45.79.80.100:data/ensembl82/fasta/ ~/data/ensembl82/fasta
```

### gff3


```bash
mkdir -p ~/data/ensembl82/gff3
cd ~/data/ensembl82/gff3

# All plants
wget -m ftp://ftp.ensemblgenomes.org/pub/plants/release-29/gff3 .

mv ftp.ensemblgenomes.org/pub/plants/release-29/gff3/* .

find . -name ".listing" | xargs rm
```

## Build local databases

Use `build_ensembl.pl`.

```bash
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --checksum --ensembl ~/data/ensembl82/mysql/homo_sapiens_core_82_37

perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --initdb --db human_82 --ensembl ~/data/ensembl82/mysql/homo_sapiens_core_82_37
```

Or use `ensembl_batch.pl`, see [this](README.md).
