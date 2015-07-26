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

## Process plastid genomes

Work flow.

```text
id ---> lineage ---> filtering ---> naming ---> strain_info.pl   ---> strain_bz.pl
                                      |                                 ^
                                      |-------> batch_get_seq.pl -------|
```

### Scrap id and acc from NCBI

Open browser and visit [NCBI plastid page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=plastid).
Save page to a local file, html only. In this case, it's `doc/green_plants_plastid_150725.html`.

Viridiplantae (33090) plastid genomes - 680 records
    Chlorophyta (3041)  [46]
    Streptophyta (35493)  [634]

From now on, our cwd is `~/data/organelle/plastid.new`.

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions.

Got **678** accessions.

```bash
mkdir -p ~/data/organelle/plastid.new
cd ~/data/organelle/plastid.new

# id,acc
# 996148,NC_017006
perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl ~/Scripts/withncbi/doc/green_plants_plastid_150725.html > plastid_id_seq.csv

# 678
cat plastid_id_seq.csv | grep -v "^#" | wc -l
```

### Add linage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `plastid.CHECKME.csv`.

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

# Koliella corcontica (a green algae) was grouped to Streptophyta.
# Don't fix it here.
```

### Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

We got 71 genera,  322 species and **352** accessions.

```text
678 ---------> 675 ---------> 352
        NA           genus
```

```bash
# filter out accessions without linage information
cat plastid.CHECKME.csv \
    | perl -nl -a -F"," -e \
    '/^#/ and next; ($F[3] eq q{NA} or $F[4] eq q{NA} ) and next; print' \
    > plastid.tmp

wc -l plastid.tmp

# valid genera
cat plastid.tmp \
    | perl -nl -a -F"," -e \
    '$seen{$F[4]}++; END {for $k (sort keys %seen) {printf qq{,%s,\n}, $k if $seen{$k} > 1}}' \
    > genus.tmp

# intersect between two files
grep -F -f genus.tmp plastid.tmp > plastid.genus.tmp
wc -l plastid.genus.tmp

# count every ranks
#   82 genus.txt
#  345 species.txt
cut -d',' -f 5 plastid.genus.tmp | sort | uniq > genus.txt
cut -d',' -f 4 plastid.genus.tmp | sort | uniq > species.txt
wc -l genus.txt species.txt

# results produced in this step
head -n 1 plastid.CHECKME.csv > plastid.FILTERED.csv
cat plastid.genus.tmp >> plastid.FILTERED.csv

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
# Fagopyrum esculentum subsp. ancestrale
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
```

Create abbreviations.

```bash
# # use Text::Abbrev
# cat genus.txt \
#     | perl -MText::Abbrev -nl -e \
#     'push @list, $_; END{%hash = abbrev @list; @ks = sort keys %hash; for $i (reverse(0 .. $#ks)) {if (index($ks[$i], $ks[$i - 1]) != 0) { print $ks[$i], q{,}, $hash{$ks[$i]}; } } }' \
#     | perl -e 'print reverse <>'
#
# # use MyUtil::abbr_most
# cat genus.txt \
#     | perl -I ~/Scripts/withncbi/lib -MMyUtil -MYAML -nl -e \
#     'push @list, $_; END{$hashref = MyUtil::abbr_most( \@list); print Dump $hashref }'


echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr' > plastid.ABBR.csv
cat plastid.FILTERED.csv \
    | grep -v '^#' \
    | perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 \
    >> plastid.ABBR.csv
```

### Download sequences and regenerate lineage information.

We don't rename sequences here, so the file has three columns

```bash
mkdir -p ~/data/organelle/plastid.new
cd ~/data/organelle/plastid.new

echo "#strain_name,accession,strain_taxon_id" > plastid_name_acc_id.csv
cat plastid.ABBR.csv \
    | grep -v '^#' \
    | perl -nl -a -F"," -e 'print qq{$F[9],$F[1],$F[0]}' \
    | sort \
    >> plastid_name_acc_id.csv

# some warnings fro bioperl, normally just ignore them
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl -f plastid_name_acc_id.csv -p 2>&1 | tee plastid_seq.log

# rsync --progress -av wangq@139.162.23.84:/home/wangq/data/organelle/ ~/data/organelle/

# count downloaded sequences
find . -name "*.fasta" | wc -l
```

Regenerate `plastid_ncbi.csv` with abbr names.

```bash
cat plastid.ABBR.csv \
    | grep -v '^#' \
    | perl -nl -a -F"," -e 'print qq{$F[0],$F[9]}' \
    | perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file plastid_ncbi.csv
```
