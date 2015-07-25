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

### Scrap id and acc from NCBI

Open browser and visit [NCBI plastid page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=plastid).
Save page to a local file, html only. In this case, it's `doc/green_plants_plastid_150725.html`.

```bash
mkdir -p ~/data/organelle/plastid.new
cd ~/data/organelle/plastid.new

# extract Taxonomy ids and genbank accessions.
# id,acc
# 996148,NC_017006
perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl ~/Scripts/withncbi/doc/green_plants_plastid_150725.html > plastid_id_seq.csv

# 678
cat plastid_id_seq.csv | grep -v "^#" | wc -l
```

### Add linage information

Got **678** accessions, give them better shapes for manually checking and automatic filtering.

```bash
# generate a .csv file for manually checking
echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum' > plastid.CHECKME.csv
cat plastid_id_seq.csv \
    | grep -v strain_taxon_id \
    | perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," \
    | perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank species \
    | perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank genus \
    | perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank family \
    | perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank order \
    | perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank class \
    | perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank phylum \
    >> plastid.CHECKME.csv

# correct 'Chlorella' mirabilis  
# darwin (bsd) need "" for -i
sed -i "" "s/\'//g" plastid.CHECKME.csv
```

### Filtering based on valid families and genera

Family has 3 or more strains and Genus has 2 or more.

We got 31 families, 71 genera,  322 species and **329** accessions.

```
        NA           family         genus
678 ---------> 649 ---------> 543 ---------> 329
```

```bash
# filter out accessions without linage information
cat plastid.CHECKME.csv \
    | perl -nl -a -F"," -e \
    '/^#/ and next; $F[5] eq q{NA} and next; $F[4] eq q{NA} and next; $F[3] eq q{NA} and next; print' \
    > plastid.tmp

# 649
wc -l plastid.tmp

# valid families
cat plastid.tmp \
    | perl -nl -a -F"," -e \
    '$seen{$F[5]}++; END {for $k (sort keys %seen) {printf qq{,%s,\n}, $k if $seen{$k} > 2}}' \
    > family.tmp

# 543
grep -F -f family.tmp plastid.tmp > plastid.family.tmp
wc -l plastid.family.tmp

# valid genera
cat plastid.family.tmp \
    | perl -nl -a -F"," -e \
    '$seen{$F[4]}++; END {for $k (sort keys %seen) {printf qq{,%s,\n}, $k if $seen{$k} > 1}}' \
    > genus.tmp

# 329
grep -F -f genus.tmp plastid.family.tmp > plastid.familiy.genus.tmp
wc -l plastid.familiy.genus.tmp

# count every ranks
#   31 family.txt
#   71 genus.txt
#  322 species.txt
cut -d',' -f 6 plastid.familiy.genus.tmp | sort | uniq > family.txt
cut -d',' -f 5 plastid.familiy.genus.tmp | sort | uniq > genus.txt
cut -d',' -f 4 plastid.familiy.genus.tmp | sort | uniq > species.txt
wc -l family.txt genus.txt species.txt

# results produced in this step
head -n 1 plastid.CHECKME.csv > plastid.FILTERED.csv
cat plastid.familiy.genus.tmp >> plastid.FILTERED.csv

# clean
rm *.tmp
```

### Find a way to name these.

 Seems it's OK to use species as names.

```bash
# sub-species
cat plastid.FILTERED.csv \
    | perl -nl -a -F"," -e \
    '/^#/i and next; $seen{$F[3]}++; END {for $k (keys %seen){printf qq{%s,%d\n}, $k, $seen{$k} if $seen{$k} > 1}};' \
    | sort

# Fragaria vesca,2
# Magnolia officinalis,2
# Olea europaea,4
# Oryza sativa,2
# Saccharum hybrid cultivar,2

# strain name not equal to species
cat plastid.FILTERED.csv \
    | grep -v '^#' \
    | perl -nl -a -F"," -e '$F[2] ne $F[3] and print $F[2]' \
    | sort

# Brassica rapa subsp. pekinensis
# Cucumis melo subsp. melo
# Fragaria vesca subsp. bracteata
# Fragaria vesca subsp. vesca
# Hordeum vulgare subsp. vulgare
# Magnolia officinalis subsp. biloba
# Olea europaea subsp. cuspidata
# Olea europaea subsp. europaea
# Olea europaea subsp. maroccana
# Oryza sativa Indica Group
# Oryza sativa Japonica Group
# Phalaenopsis aphrodite subsp. formosana
# Saccharum hybrid cultivar NCo 310
# Saccharum hybrid cultivar SP80-3280

# create abbreviations

# # use Text::Abbrev
# cat genus.txt \
#     | perl -MText::Abbrev -nl -e \
#     'push @list, $_; END{%hash = abbrev @list; @ks = sort keys %hash; for $i (reverse(0 .. $#ks)) {if (index($ks[$i], $ks[$i - 1]) != 0) { print $ks[$i], q{,}, $hash{$ks[$i]}; } } }' \
#     | perl -e 'print reverse <>' \
#     > genus.abbr.txt
#
# # use MyUtil::abbr_most
# cat genus.txt \
#     | perl -I ~/Scripts/withncbi/lib -MMyUtil -MYAML -nl -e \
#     'push @list, $_; END{$hashref = MyUtil::abbr_most( \@list); print Dump $hashref }'

cat plastid.FILTERED.csv \
    | grep -v '^#' \
    | perl ~/Scripts/withncbi/taxon/auto_name.pl -c "3,4,5" -s ","
```

```bash
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl plastid_id_seq.csv  2>&1 | tee plastid_seq.log

# count downloaded sequences
find . -name "*.fasta" | wc -l

```
