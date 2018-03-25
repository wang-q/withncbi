# Process plant mitochondrion genomes

[TOC levels=1-3]: # " "
- [Process plant mitochondrion genomes](#process-plant-mitochondrion-genomes)
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
- [Add lineage information](#add-lineage-information)
    - [Can't get clear taxon information](#cant-get-clear-taxon-information)
- [Filtering based on valid families and genera](#filtering-based-on-valid-families-and-genera)
- [Find a way to name these.](#find-a-way-to-name-these)
- [Download sequences and regenerate lineage information.](#download-sequences-and-regenerate-lineage-information)
- [Create alignment plans](#create-alignment-plans)


# Scrap id and acc from NCBI

Open browser and visit
[NCBI mitochondrion page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=organelle).
Save page to a local file, html only. In this case, it's
`doc/green_plants_mitochondrion_180325.html`.

All [Eukaryota](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=organelle),
`doc/eukaryota_mitochondrion_180325.html`.

```text
Eukaryota (2759)                8455
    Viridiplantae (33090)       211
        Chlorophyta (3041)      51
        Streptophyta (35493)    160
```

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions from all history
pages.

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

# rice mitochondria have plasmids
sed -i".bak" "s/,NC_001751$/,NC_011033/" plant_mitochondrion_id_seq.csv # japonica
sed -i".bak" "s/,NC_001776$/,NC_007886/" plant_mitochondrion_id_seq.csv # indica

```

# Add lineage information

Give ids better shapes for manually checking and automatic filtering.

*Update `~/data/NCBI/taxdmp` before running `id_project_to.pl`*.

If you sure, you can add or delete lines and contents in `mitochondrion.CHECKME.csv`.

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


Split Streptophyta according to http://www.theplantlist.org/

```bash
cd ~/data/organelle/mitochondrion_summary

# Angiosperms
perl -Mojo -e '
    g(q{http://www.theplantlist.org/browse/A/})->dom
    ->find(q{li > a > i[class=family]})
    ->each( sub { print shift->text . "\n" } );
    ' > Angiosperms.tmp
echo Aceraceae >> Angiosperms.tmp
echo Asteraceae >> Angiosperms.tmp
echo Campynemataceae >> Angiosperms.tmp
echo Chenopodiaceae >> Angiosperms.tmp
echo Fabaceae >> Angiosperms.tmp
echo Viscaceae >> Angiosperms.tmp

cat Angiosperms.tmp |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Angiosperms/g
        '\'' mitochondrion.CHECKME.csv
    '

# Gymnosperms
perl -Mojo -e '
    g(q{http://www.theplantlist.org/browse/G/})->dom
    ->find(q{li > a > i[class=family]})
    ->each( sub { print shift->text . "\n" } );
    ' |
    (echo Sciadopityaceae && cat) |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Gymnosperms/g
        '\'' mitochondrion.CHECKME.csv
    '

# Pteridophytes
perl -Mojo -e '
    g(q{http://www.theplantlist.org/browse/P/})->dom
    ->find(q{li > a > i[class=family]})
    ->each( sub { print shift->text . "\n" } );
    ' |
    (echo Lygodiaceae && cat) |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Pteridophytes/g
        '\'' mitochondrion.CHECKME.csv
    '

# Bryophytes
perl -Mojo -e '
    g(q{http://www.theplantlist.org/browse/B/})->dom
    ->find(q{li > a > i[class=family]})
    ->each( sub { print shift->text . "\n" } );
    ' |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Bryophytes/g
        '\'' mitochondrion.CHECKME.csv
    '

rm *.tmp *.bak

```

## Can't get clear taxon information

FIXME

# Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

```text
217 ---------> 209 ---------> 89 ---------> 111
        NA           genus         family
```


```bash
mkdir -p ~/data/organelle/mitochondrion_summary
cd ~/data/organelle/mitochondrion_summary

# filter out accessions without linage information (strain, species, genus and family)
cat mitochondrion.CHECKME.csv |
    perl -nla -F"," -e '
        /^#/ and next;
        ($F[2] eq q{NA} or $F[3] eq q{NA} or $F[4] eq q{NA} or $F[5] eq q{NA} ) and next;
        print
    ' \
    > mitochondrion.tmp

# 209
wc -l mitochondrion.tmp

#----------------------------#
# Genus
#----------------------------#
# valid genera
cat mitochondrion.tmp |
    perl -nla -F"," -e '
        $seen{$F[4]}++; 
        END {
            for $k (sort keys %seen) {
                printf qq{,%s,\n}, $k if $seen{$k} > 1
            }
        }
    ' \
    > genus.tmp

# intersect between two files
grep -F -f genus.tmp mitochondrion.tmp > mitochondrion.genus.tmp

# 89
wc -l mitochondrion.genus.tmp

#----------------------------#
# Family
#----------------------------#
# get some genera back as candidates for outgroup
cat mitochondrion.genus.tmp |
    perl -nla -F"," -e 'printf qq{,$F[5],\n}' \
    > family.tmp

# intersect between two files
grep -F -f family.tmp mitochondrion.tmp > mitochondrion.family.tmp

# 111
wc -l mitochondrion.family.tmp

#----------------------------#
# results produced in this step
#----------------------------#
head -n 1 mitochondrion.CHECKME.csv > mitochondrion.DOWNLOAD.csv
cat mitochondrion.family.tmp >> mitochondrion.DOWNLOAD.csv

# clean
rm *.tmp *.bak
```

# Find a way to name these.

Seems it's OK to use species as names.

```bash
# sub-species
cat mitochondrion.DOWNLOAD.csv |
    perl -nl -a -F"," -e '
        /^#/i and next; 
        $seen{$F[3]}++; 
        END {
            for $k (keys %seen){printf qq{%s,%d\n}, $k, $seen{$k} if $seen{$k} > 1}
        };
    ' |
    sort

#Beta vulgaris,2
#Oryza sativa,2
#Zea mays,3

# strain name not equal to species
cat mitochondrion.DOWNLOAD.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e '$F[2] ne $F[3] and print $F[2]' |
    sort

#Beta vulgaris subsp. maritima
#Beta vulgaris subsp. vulgaris
#Brassica rapa subsp. oleifera
#Calypogeia fissa subsp. neogaea
#Oryza sativa Indica Group
#Oryza sativa Japonica Group
#Zea mays subsp. mays
#Zea mays subsp. parviglumis

```

Create abbreviations.

```bash
cd ~/data/organelle/mitochondrion_summary

echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr' > mitochondrion.ABBR.csv
cat mitochondrion.DOWNLOAD.csv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 |
    sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> mitochondrion.ABBR.csv

```

# Download sequences and regenerate lineage information.

```bash
mkdir -p ~/data/organelle/mitochondrion_genomes
cd ~/data/organelle/mitochondrion_genomes

echo "#strain_name,accession,strain_taxon_id" > mitochondrion_name_acc_id.csv
cat ../mitochondrion_summary/mitochondrion.ABBR.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e 'print qq{$F[9],$F[1],$F[0]}' |
    sort \
    >> mitochondrion_name_acc_id.csv

# local, Runtime 10 seconds.
# with --entrez, Runtime 7 minutes and 23 seconds.
# And which-can't-find is still which-can't-find.
cat ../mitochondrion_summary/mitochondrion.ABBR.csv |
    grep -v '^#' |
    perl -nla -F"," -e 'print qq{$F[0],$F[9]}' |
    uniq |
    perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file mitochondrion_ncbi.csv

# some warnings from bioperl, just ignore them
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -f mitochondrion_name_acc_id.csv \
    -p 2>&1 |
    tee mitochondrion_seq.log

# count downloaded sequences
find . -name "*.fasta" | wc -l

```


# Create alignment plans

We got **111** accessions.

Numbers for higher ranks are: 15 orders, 17 families, 27 genera and 85 species.

```bash
cd ~/data/organelle/mitochondrion_summary

# valid genera
cat mitochondrion.ABBR.csv |
    grep -v "^#" |
    perl -nl -a -F"," -e '
        $seen{$F[4]}++; 
        END {
            for $k (sort keys %seen) {
                printf qq{,%s,\n}, $k if $seen{$k} > 1
            }
        }
    ' \
    > genus.tmp

# intersect between two files
grep -F -f genus.tmp mitochondrion.ABBR.csv > mitochondrion.GENUS.csv

# 89
wc -l mitochondrion.GENUS.csv

#   count every ranks
#      15 order.list.tmp
#      17 family.list.tmp
#      27 genus.list.tmp
#      85 species.list.tmp
cut -d',' -f 4 mitochondrion.GENUS.csv | sort | uniq > species.list.tmp
cut -d',' -f 5 mitochondrion.GENUS.csv | sort | uniq > genus.list.tmp
cut -d',' -f 6 mitochondrion.GENUS.csv | sort | uniq > family.list.tmp
cut -d',' -f 7 mitochondrion.GENUS.csv | sort | uniq > order.list.tmp
wc -l order.list.tmp family.list.tmp genus.list.tmp species.list.tmp

# create again with headers
grep -F -f genus.tmp mitochondrion.ABBR.csv > mitochondrion.GENUS.tmp

# sort by multiply columns, phylum, order, family, abbr
head -n 1 mitochondrion.ABBR.csv > mitochondrion.GENUS.csv
cat mitochondrion.GENUS.tmp \
    | sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> mitochondrion.GENUS.csv

# clean
rm *.tmp *.bak
```

