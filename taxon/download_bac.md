# Bulk downloading bacteria genomes from NCBI

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
