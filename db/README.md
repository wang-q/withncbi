# `db/`

Turn NCBI genome reports and assembly reports into query-able MySQL databases.

Also, taxonomy information is added to all items.

Downloading date: 2020-12-11

[TOC levels=1-3]: # ""

- [`db/`](#db)
- [Get data from NCBI](#get-data-from-ncbi)
- [Databases](#databases)
  - [Genome reports](#genome-reports)
  - [Assembly reports](#assembly-reports)
- [Old Bacteria genomes](#old-bacteria-genomes)


# Get data from NCBI

Paths in the NCBI ftp site:

* <ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS>
* <ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS>

Local paths listed in `config.ini`.

NCBI also provides other download methods, including `rsync` and `aspera`.

I use the following command lines on a Linux box. For mac, aspera's path is different.

```shell script
# gr
rsync -avP ftp.ncbi.nlm.nih.gov::genomes/GENOME_REPORTS/ \
    ~/data/NCBI/genomes/GENOME_REPORTS/

#~/.aspera/connect/bin/ascp \
#    -TQ -k1 -p -v \
#    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
#    anonftp@ftp-private.ncbi.nlm.nih.gov:/genomes/ASSEMBLY_REPORTS \
#   ~/data/NCBI/genomes/

# ar
# there're hidden useless directories.
rsync -avP ftp.ncbi.nlm.nih.gov::genomes/ASSEMBLY_REPORTS/ \
    --exclude=".tmp" \
    --exclude=".old" \
    ~/data/NCBI/genomes/ASSEMBLY_REPORTS/

```

NCBI bioproject and taxonomy are also needed.

```shell script
# bioproject
rsync -avP ftp.ncbi.nlm.nih.gov::bioproject/ \
    --exclude="*.xml" \
    ~/data/NCBI/bioproject/

# taxonomy
rsync -avP ftp.ncbi.nlm.nih.gov::pub/taxonomy/ \
    --exclude=".tmp" \
    --exclude=".old" \
    --exclude="*.Z" \
    --exclude="taxdump_archive" \
    --exclude="new_taxdump" \
    --exclude="accession2taxid" \
    --exclude="gi_taxid_*" \
    ~/data/NCBI/taxonomy/

rm -fr ~/data/NCBI/taxdmp
mkdir -p ~/data/NCBI/taxdmp
tar xvfz ~/data/NCBI/taxonomy/taxdump.tar.gz -C ~/data/NCBI/taxdmp

```

# Databases

We will create 4 MySQL databases:

* gr_prok: genome reports for prokaryotes;
* gr_euk: genome reports for eukaryotes;
* ar_refseq: assembly reports for RefSeq;
* ar_genbank: assembly reports for GenBank.

Also, generate some useful excel workbooks.

## Genome reports

```shell script
cd ~/Scripts/withncbi/db

# raw genome reports to .csv files
perl gr_strains.pl -o prok_strains.csv
perl gr_strains.pl --euk -o euk_strains.csv

# load .csv to MySQL
perl gr_db.pl --db gr_prok --file prok_strains.csv
perl gr_db.pl --db gr_euk --file euk_strains.csv

# generate .xlsx
perl gr_overview.pl --db gr_prok
perl gr_overview.pl --db gr_euk

```

## Assembly reports

```shell script
cat ~/data/NCBI/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt \
    ~/data/NCBI/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt |
    grep -v "^#" |
    tsv-select -f 12 | # assembly_level
    sort |
    uniq -c
#13066 Chromosome
#90279 Complete Genome
#785524 Contig
#215038 Scaffold

```

```shell script
cd ~/Scripts/withncbi/db

perl ar_strains.pl -o ar_strains.csv
perl ar_strains.pl --genbank -o ar_strains_genbank.csv

wc -l *strains*.csv
#   37480 ar_strains.csv
#   46042 ar_strains_genbank.csv
#    6761 euk_strains.csv
#   37471 prok_strains.csv

perl ar_db.pl --db ar_refseq --file ar_strains.csv
perl ar_db.pl --db ar_genbank --file ar_strains_genbank.csv

perl ar_overview.pl --db ar_refseq
perl ar_overview.pl --db ar_genbank

#cp -f *.xlsx ../doc
#rm *.xlsx *.csv

```


# Old Bacteria genomes

On 02 December 2015, these directories were moved to
`ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/`.

```shell script
rsync -av -P ftp.ncbi.nlm.nih.gov::genomes/archive/old_refseq/Bacteria/ \
    --exclude="all.*" \
    --exclude=".tmp" \
    --exclude=".old" \
    ~/data/NCBI/genomes/Bacteria/

rsync -av -P ftp.ncbi.nlm.nih.gov::genomes/archive/old_refseq/Bacteria_DRAFT/ \
    --exclude=".tmp" \
    --exclude=".old" \
    ~/data/NCBI/genomes/Bacteria_DRAFT/

rsync -av -P ftp.ncbi.nlm.nih.gov::genomes/archive/old_genbank/Bacteria/ \
    --exclude=".tmp" \
    --exclude=".old" \
    ~/data/NCBI/genbank/genomes/Bacteria/

rsync -av -P ftp.ncbi.nlm.nih.gov::genomes/archive/old_genbank/Bacteria_DRAFT/ \
    --exclude=".tmp" \
    --exclude=".old" \
    ~/data/NCBI/genbank/genomes/Bacteria_DRAFT/

```

Newer genomes list in genomes/refseq/bacteria are just symlinks to genomes/all/*.

So local mirrors are no longer needed.
