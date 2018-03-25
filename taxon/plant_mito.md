# Process plant mitochondrion genomes

[TOC]: # " "
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
- [Add lineage information](#add-lineage-information)


## Scrap id and acc from NCBI

Open browser and visit
[NCBI plastid page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=organelle).
Save page to a local file, html only. In this case, it's
`doc/green_plants_mitochondrion_180325.html`.

All
[Eukaryota](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=organelle),
`doc/eukaryota_mitochondrion_180325.html`.

```text
Eukaryota (2759)                8455
    Viridiplantae (33090)       211
        Chlorophyta (3041)      51
        Streptophyta (35493)    160
```

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank
accessions from all history pages.

Got **217** accessions.

```bash
mkdir -p ~/data/organelle/mitochondrion_genomes
cd ~/data/organelle/mitochondrion_genomes

rm webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/eukaryota_mitochondrion_180325.html \
    >> webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/green_plants_mitochondrion_180325.html \
    >> webpage_id_seq.csv    

```

Use `taxon/gb_taxon_locus.pl` to extract information from refseq genbank files.

```bash
cd ~/data/organelle/mitochondrion_genomes

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz

gzip -c -d mitochondrion.*.genomic.gbff.gz > mitochondrion.genomic.gbff

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl mitochondrion.genomic.gbff > refseq_id_seq.csv

rm mitochondrion.genomic.gbff

# 8566
cat refseq_id_seq.csv | grep -v "^#" | wc -l

# combine
cat webpage_id_seq.csv refseq_id_seq.csv |
    sort -u -t, -k1,1 \
    > mitochondrion_id_seq.csv

# 8535
cat mitochondrion_id_seq.csv | grep -v "^#" | wc -l

```

Restrict taxonomy ids to green plants with `taxon/id_restrict.pl`.

```bash
cd ~/data/organelle/mitochondrion_genomes

echo '#strain_taxon_id,accession' > plant_mitochondrion_id_seq.csv
cat mitochondrion_id_seq.csv |
    grep -v "^#" |
    perl ~/Scripts/withncbi/taxon/id_restrict.pl -s "," -a 33090 \
    >> plant_mitochondrion_id_seq.csv

# 217
cat plant_mitochondrion_id_seq.csv | grep -v "^#" | wc -l

```

## Add lineage information

Give ids better shapes for manually checking and automatic filtering.

*Update `~/data/NCBI/taxdmp` before running `id_project_to.pl`*.

If you sure, you can add or delete lines and contents in
`mitochondrion.CHECKME.csv`.

```bash
mkdir -p ~/data/organelle/mitochondrion_summary
cd ~/data/organelle/mitochondrion_summary

# generate a .csv file for manually checking
echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum' > mitochondrion.CHECKME.csv
cat ../mitochondrion_genomes/plant_mitochondrion_id_seq.csv |
    grep -v "^#" |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank species |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank genus |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank family |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank order |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank class |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank phylum |
    sort -t',' -k9,9 -k8,8 -k7,7 -k6,6 -k5,5 \
    >> mitochondrion.CHECKME.csv

```

Manually correct lineages. # FIXME

