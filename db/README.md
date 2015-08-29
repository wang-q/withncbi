# `db/`

Turn NCBI genome reports and assembly reports into query-able MySQL databases.

Also, taxonomy information are added to all items.

## Get data from NCBI

Download paths from NCBI ftp:

* ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS
* ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS

Local paths listed in `config.ini`.  

NCBI also provides other download methods including rsync and aspera.

I use the following command lines on a linux box. For mac, aspera's path is different.

```bash
# gr
rsync --progress -av ftp.ncbi.nlm.nih.gov::genomes/GENOME_REPORTS/ \
    ~/data/NCBI/genomes/GENOME_REPORTS/

# ar is huge
~/.aspera/connect/bin/ascp \
    -TQ -k1 -p -v \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    anonftp@ftp-private.ncbi.nlm.nih.gov:/genomes/ASSEMBLY_REPORTS \
   ~/data/NCBI/genomes/

# there're hidden useless directories.
rsync --progress -av ftp.ncbi.nlm.nih.gov::genomes/ASSEMBLY_REPORTS/All/ \
    ~/data/NCBI/genomes/ASSEMBLY_REPORTS/All/
```

NCBI bioproject and taxonomy is also needed.

```bash
# bioproject
rsync --progress -av ftp.ncbi.nlm.nih.gov::bioproject/ \
    ~/data/NCBI/bioproject/

# taxonomy
rsync --progress -av ftp.ncbi.nlm.nih.gov::pub/taxonomy/ \
    ~/data/NCBI/taxonomy/
```

Bacteria genomes.

```bash
rsync -av -P ftp.ncbi.nlm.nih.gov::genomes/Bacteria/ \
    ~/data/NCBI/genomes/Bacteria/ \
    --exclude="all.*"

rsync -av -P ftp.ncbi.nlm.nih.gov::genomes/Bacteria_DRAFT/ \
    ~/data/NCBI/genomes/Bacteria_DRAFT/

rsync -av -P ftp.ncbi.nlm.nih.gov::genbank/genomes/Bacteria/ \
    ~/data/NCBI/genbank/genomes/Bacteria/

rsync -av -P ftp.ncbi.nlm.nih.gov::genbank/genomes/Bacteria_DRAFT/ \
    ~/data/NCBI/genbank/genomes/Bacteria_DRAFT/

```

## Databases

We create 4 MySQL databases:

    * gr_prok: genome reports for prokaryotes;
    * gr_euk: genome reports for eukaryotes;
    * ar_refseq: assembly reports for RefSeq;
    * ar_genbank: assembly reports for GenBank.

Also generate some useful excel workbooks.

### Genome reports

```bash
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

### Assembly reprots

```bash
cd ~/Scripts/withncbi/db

perl ar_strains.pl -o ar_strains.csv
perl ar_strains.pl --genbank -o ar_strains_genbank.csv

perl ar_db.pl --db ar_refseq --file ar_strains.csv
perl ar_db.pl --db ar_genbank --file ar_strains_genbank.csv

perl ar_overview.pl --db ar_refseq
perl ar_overview.pl --db ar_genbank

cp -f *.xlsx ../db
rm *.xlsx *.csv
```
