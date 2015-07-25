## Bulk downloading bacteria genomes from NCBI

On a linux box.

```bash
# genome Bacteria
rsync --progress -av ftp.ncbi.nlm.nih.gov::genomes/Bacteria/ \
    ~/data/NCBI/genomes/Bacteria/ --exclude="all.*"

# genomes Bacteria_DRAFT
~/.aspera/connect/bin/ascp \
    -TQ -k1 -p -v \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    anonftp@ftp-private.ncbi.nlm.nih.gov:/genomes/Bacteria_DRAFT \
   ~/data/NCBI/genomes/

rsync --progress -av ftp.ncbi.nlm.nih.gov::genomes/Bacteria_DRAFT/ \
    ~/data/NCBI/genomes/Bacteria_DRAFT/
```

NCBI will abandon genbank Bacteria and genbank Bacteria_DRAFT soon. And updates have stopped in late
2014,

ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/NOTICE_OF_CHANGE.txt

```bash
# genbank Bacteria
~/.aspera/connect/bin/ascp \
    -TQ -k1 -p -v \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    anonftp@ftp-private.ncbi.nlm.nih.gov:/genbank/genomes/Bacteria \
   ~/data/NCBI/genbank/genomes/

rsync --progress -av ftp.ncbi.nlm.nih.gov::genbank/genomes/Bacteria/ \
    ~/data/NCBI/genbank/genomes/Bacteria/

# genbank Bacteria_DRAFT
~/.aspera/connect/bin/ascp \
    -TQ -k1 -p -v \
    -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    anonftp@ftp-private.ncbi.nlm.nih.gov:/genbank/genomes/Bacteria_DRAFT \
   ~/data/NCBI/genbank/genomes/

rsync --progress -av ftp.ncbi.nlm.nih.gov::genbank/genomes/Bacteria_DRAFT/ \
    ~/data/NCBI/genbank/genomes/Bacteria_DRAFT/
```

## Download plastid genomes

Open browser and visit [NCBI plastid page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=plastid).
Save page to local file, html only. In this time, it's `doc/Eukaryota_plastid_150725.html`.

```bash
mkdir -p ~/data/organelle/plastid.new
cd ~/data/organelle/plastid.new

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl ~/Scripts/withncbi/doc/Eukaryota_plastid_150725.html > plastid_id_seq.csv
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl plastid_id_seq.csv  2>&1 | tee plastid_seq.log

# count downloaded sequences
find . -name "*.fasta" | wc -l

```
