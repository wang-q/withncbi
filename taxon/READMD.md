# Process plastid genomes

The following command lines are about how I processed the plastid genomes of green plants.
Many tools of `taxon/` are used here, which makes a good example for users.

## Work flow.

```text
id ---> lineage ---> filtering ---> naming ---> strain_info.pl   ---> strain_bz.pl
                                      |                                 ^
                                      |-------> batch_get_seq.pl -------|
```

I'm sure there are no commas in names. So for convenient, don't use Text::CSV_XS.

## Scrap id and acc from NCBI

Open browser and visit [NCBI plastid page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=plastid).
Save page to a local file, html only. In this case, it's `doc/green_plants_plastid_150725.html`.

```text
Viridiplantae (33090) plastid genomes - 680 records
    Chlorophyta (3041)  [46]
    Streptophyta (35493)  [634]
```

From now on, our cwd is `~/data/organelle/plastid_genomes`.

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions.

Got **678** accessions.

```bash
mkdir -p ~/data/organelle/plastid_genomes
cd ~/data/organelle/plastid_genomes

# id,acc
# 996148,NC_017006
perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl ~/Scripts/withncbi/doc/green_plants_plastid_150725.html > plastid_id_seq.csv

# 678
cat plastid_id_seq.csv | grep -v "^#" | wc -l
```

## Add linage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `plastid.CHECKME.csv`.

```bash
# generate a .csv file for manually checking
echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum' > plastid.CHECKME.csv
cat plastid_id_seq.csv \
    | grep -v "^#" \
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
sed -i ".bak" "s/\'//g" plastid.CHECKME.csv

sed -i ".bak" "s/Chlorella,NA,NA/Chlorella,Chlorellaceae,Chlorellales/" plastid.CHECKME.csv

# Koliella corcontica (a green algae) was grouped to Streptophyta.
sed -i ".bak" "s/Klebsormidiophyceae,Streptophyta/Klebsormidiophyceae,Chlorophyta/" plastid.CHECKME.csv

sed -i ".bak" "s/Koliella,NA/Koliella,Klebsormidiaceae/" plastid.CHECKME.csv

sed -i ".bak" "s/Nephroselmis,NA,NA/Nephroselmis,Nephroselmidaceae,Nephroselmidales/" plastid.CHECKME.csv

# Chrysanthemum x morifolium and Pelargonium x hortorum are also weird, but they can be googled.
```

## Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

```text
678 ---------> 675 ---------> 352 ---------> 526
        NA           genus          family
```

```bash
# filter out accessions without linage information
cat plastid.CHECKME.csv \
    | perl -nl -a -F"," -e \
    '/^#/ and next; ($F[3] eq q{NA} or $F[4] eq q{NA} ) and next; print' \
    > plastid.tmp

# 675
wc -l plastid.tmp

#----------------------------#
# Genus
#----------------------------#
# valid genera
cat plastid.tmp \
    | perl -nl -a -F"," -e \
    '$seen{$F[4]}++; END {for $k (sort keys %seen) {printf qq{,%s,\n}, $k if $seen{$k} > 1}}' \
    > genus.tmp

# intersect between two files
grep -F -f genus.tmp plastid.tmp > plastid.genus.tmp

# 352
wc -l plastid.genus.tmp

#----------------------------#
# Family
#----------------------------#
# get some genera back as candidates for outgroup
cat plastid.genus.tmp \
    | perl -nl -a -F"," -e 'printf qq{,$F[5],\n}' \
    > family.tmp

# intersect between two files
grep -F -f family.tmp plastid.tmp > plastid.family.tmp

# 526
wc -l plastid.family.tmp

# results produced in this step
head -n 1 plastid.CHECKME.csv > plastid.DOWNLOAD.csv
cat plastid.family.tmp >> plastid.DOWNLOAD.csv

# clean
rm *.tmp *.bak
```

## Find a way to name these.

Seems it's OK to use species as names.

```bash
# sub-species
cat plastid.DOWNLOAD.csv \
    | perl -nl -a -F"," -e \
    '/^#/i and next; $seen{$F[3]}++; END {for $k (keys %seen){printf qq{%s,%d\n}, $k, $seen{$k} if $seen{$k} > 1}};' \
    | sort

# Fragaria vesca,2
# Magnolia officinalis,2
# Olea europaea,4
# Oryza sativa,2
# Saccharum hybrid cultivar,2

# strain name not equal to species
cat plastid.DOWNLOAD.csv \
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
cat plastid.DOWNLOAD.csv \
    | grep -v '^#' \
    | perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 \
    >> plastid.ABBR.csv
```

## Download sequences and regenerate lineage information.

We don't rename sequences here, so the file has three columns

```bash
mkdir -p ~/data/organelle/plastid_genomes
cd ~/data/organelle/plastid_genomes

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

Create `plastid_ncbi.csv` with abbr names.

```bash
# local, Runtime 10 seconds.
# with --entrez, Runtime 7 minutes and 23 seconds.
# And can't find is still can't find.
cat plastid.ABBR.csv \
    | grep -v '^#' \
    | perl -nl -a -F"," -e 'print qq{$F[0],$F[9]}' \
    | perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file plastid_ncbi.csv
```

## Create alignment plans

We got 43 families, 82 genera, 345 species and **352** accessions.

```bash
# valid genera
cat plastid.ABBR.csv \
    | grep -v "^#" \
    | perl -nl -a -F"," -e \
    '$seen{$F[4]}++; END {for $k (sort keys %seen) {printf qq{,%s,\n}, $k if $seen{$k} > 1}}' \
    > genus.tmp

# intersect between two files
grep -F -f genus.tmp plastid.ABBR.csv > plastid.GENUS.csv

# 352
wc -l plastid.GENUS.csv

#   count every rank
#   42 family.txt
#   82 genus.txt
#  345 species.txt
cut -d',' -f 4 plastid.GENUS.csv | sort | uniq > species.txt
cut -d',' -f 5 plastid.GENUS.csv | sort | uniq > genus.txt
cut -d',' -f 6 plastid.GENUS.csv | sort | uniq > family.txt
wc -l family.txt genus.txt species.txt

# create again with headers
grep -F -f genus.tmp plastid.ABBR.csv > plastid.GENUS.tmp

# sort by multiply columns, phylum, order, family, abbr
head -n 1 plastid.ABBR.csv > plastid.GENUS.csv
cat plastid.GENUS.tmp \
    | sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> plastid.GENUS.csv
```

Create alignments without outgroups.

```bash
cd ~/data/organelle/

echo -e "mkdir -p ~/data/organelle/plastid.working\ncd ~/data/organelle/plastid.working\n" > plastid.cmd.txt

# sed is bad for handling new lines
cat plastid_genomes/plastid.GENUS.csv \
    | grep -v "^#" \
    | perl -nl -a -F"," -e \
    'BEGIN{($g, @s) = ('');}; if ($F[4] ne $g) {if ($g) {print qq{\n# $g}; print qq{GENUS $g \\}; print qq{-q $_ \\} for @s;} $g = $F[4]; @s = ();} push @s, $F[9]; END {print qq{\n# $g}; print qq{GENUS $g \\}; print qq{-q $_ \\} for @s; print;}' \
    | perl -e '@ls = <>; $l = join q{}, @ls; $l =~ s/GENUS (\w+) \\\n\-q/GENUS \1 \\\n\    -t/gs; $l =~ s/\-q /    \-q /gs; $l =~ s/\\\n\n/\n\n/gs; print $l;' \
    | perl -e '@ls = <>; $l = join q{}, @ls; $l =~ s/GENUS /perl \~\/Scripts\/withncbi\/taxon\/strain_bz\.pl \\\nOPTIONS\\\n    --use_name \\\n    \-\-name /gs; print $l;' \
    | perl -e '@ls = <>; $l = join q{}, @ls; $l =~ s/OPTIONS/    \-\-file ~\/data\/organelle\/plastid_genomes\/plastid_ncbi.csv \\\nOPTIONS /gs; print $l;' \
    | perl -e '@ls = <>; $l = join q{}, @ls; $l =~ s/OPTIONS/    \-\-seq_dir ~\/data\/organelle\/plastid_genomes/gs; print $l;' \
    >> plastid.cmd.txt

```

## Batch running

The old prepare_run.sh

```bash
cd ~/data/organelle/

sh plastid.cmd.txt
cd plastid.working

# sh plastid.OG.cmd.txt
## cd plastid.OG

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/1_real_chr.sh ; \
    echo sh $d/2_file-rm.sh ; \
    echo sh $d/3_pair_cmd.sh ; \
    echo sh $d/4_rawphylo.sh ; \
    echo sh $d/5_multi_cmd.sh ; \
    echo sh $d/6_multi_db_only.sh ; \
    echo ; \
done  > runall.sh

sh runall.sh 2>&1 | tee plastid.final.log
# sh runall.sh 2>&1 | tee plastid.OG.final.log

# clean
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_raw" | xargs rm -fr
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | xargs rm -fr
find $PWD -mindepth 1 -maxdepth 2 -type d -name "rawphylo" | xargs rm -fr

find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy" | xargs rm
find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy.reduced" | xargs rm
```
