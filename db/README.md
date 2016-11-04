# `db/`

Turn NCBI genome reports and assembly reports into query-able MySQL
databases.

Also, taxonomy information are added to all items.

## Get data from NCBI

Download paths from NCBI ftp:

* ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS
* ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS

Local paths listed in `config.ini`.

NCBI also provides other download methods including rsync and aspera.

I use the following command lines on a linux box. For mac, aspera's path
is different.

```bash
# gr
rsync -avP ftp.ncbi.nlm.nih.gov::genomes/GENOME_REPORTS/ \
    ~/data/NCBI/genomes/GENOME_REPORTS/

# ar is huge
~/.aspera/connect/bin/ascp \
    -TQ -k1 -p -v \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    anonftp@ftp-private.ncbi.nlm.nih.gov:/genomes/ASSEMBLY_REPORTS \
   ~/data/NCBI/genomes/

# there're hidden useless directories.
rsync -avP ftp.ncbi.nlm.nih.gov::genomes/ASSEMBLY_REPORTS/ \
    --exclude=".tmp" \
    --exclude=".old" \
    ~/data/NCBI/genomes/ASSEMBLY_REPORTS/
```

NCBI bioproject and taxonomy is also needed.

```bash
# bioproject
rsync -avP ftp.ncbi.nlm.nih.gov::bioproject/ \
    ~/data/NCBI/bioproject/

# taxonomy
rsync -avP ftp.ncbi.nlm.nih.gov::pub/taxonomy/ \
    --exclude=".tmp" \
    --exclude=".old" \
    --exclude="accession2taxid" \
    ~/data/NCBI/taxonomy/
```

Old Bacteria genomes.

On 02 December 2015 these directories were moved to
`ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/`.

```bash
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

Newer genomes list in genomes/refseq/bacteria are just symlinks to
genomes/all/*.

So local mirrors are no longer needed.

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

### Assembly reports

```bash
cd ~/Scripts/withncbi/db

perl ar_strains.pl -o ar_strains.csv
perl ar_strains.pl --genbank -o ar_strains_genbank.csv

perl ar_db.pl --db ar_refseq --file ar_strains.csv
perl ar_db.pl --db ar_genbank --file ar_strains_genbank.csv

perl ar_overview.pl --db ar_refseq
perl ar_overview.pl --db ar_genbank

cp -f *.xlsx ../doc
rm *.xlsx *.csv
```

