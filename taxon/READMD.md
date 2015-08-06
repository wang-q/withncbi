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
Save page to a local file, html only. In this case, it's `doc/green_plants_plastid_150805.html`.

All [Eukaryota](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?opt=plastid&taxid=2759), `doc/eukaryota_plastid_150806.html`.

Or [this link](http://www.ncbi.nlm.nih.gov/genome/browse/?report=5).

```text
Eukaryota (2759)                901
    Viridiplantae (33090)       820
        Chlorophyta (3041)      58
        Streptophyta (35493)    762
```

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions.

Got **898** accessions.

```bash
mkdir -p ~/data/organelle/plastid_genomes
cd ~/data/organelle/plastid_genomes

# id,acc
# 996148,NC_017006
perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl ~/Scripts/withncbi/doc/eukaryota_plastid_150806.html > plastid_id_seq.csv

# 898
cat plastid_id_seq.csv | grep -v "^#" | wc -l
```

## Add linage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `plastid.CHECKME.csv`.

```bash
mkdir -p ~/data/organelle/plastid_summary
cd ~/data/organelle/plastid_summary

# generate a .csv file for manually checking
echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum' > plastid.CHECKME.csv
cat ../plastid_genomes/plastid_id_seq.csv \
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
sed -i".bak" "s/\'//g" plastid.CHECKME.csv
sed -i".bak" "s/Chlorella,NA,NA/Chlorella,Chlorellaceae,Chlorellales/" plastid.CHECKME.csv

# Koliella corcontica (a green algae) was grouped to Streptophyta.
sed -i".bak" "s/Klebsormidiophyceae,Streptophyta/Klebsormidiophyceae,Chlorophyta/" plastid.CHECKME.csv

# various missing families
# queried from http://www.algaebase.org/
sed -i".bak" "s/Aureococcus,NA/Aureococcus,Pelagomonadaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Aureoumbra,NA/Aureoumbra,Sarcinochrysidaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Bigelowiella,NA/Bigelowiella,Chlorarachniaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Choricystis,NA/Choricystis,Coccomyxaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Cryptoglena,NA,NA/Cryptoglena,Euglenaceae,Euglenales/" plastid.CHECKME.csv
sed -i".bak" "s/Dicloster,NA,NA/Dicloster,Chlorellaceae,Chlorellales/" plastid.CHECKME.csv
sed -i".bak" "s/Dictyochloropsis,NA/Dictyochloropsis,Trebouxiaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Euglenaformis,NA/Euglenaformis,Euglenaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Euglenaria,NA/Euglenaria,Euglenaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Eutreptiella,NA/Eutreptiella,Eutreptiaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Fusochloris,NA/Fusochloris,Microthamniaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Geminella,NA,NA/Geminella,Chlorellaceae,Chlorellales/" plastid.CHECKME.csv
sed -i".bak" "s/Gloeotilopsis,NA/Gloeotilopsis,Ulotrichaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Helicosporidium,NA/Helicosporidium,Chlorellaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Koliella,NA/Koliella,Klebsormidiaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Microthamnion,NA/Microthamnion,Microthamniaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Monomorphina,NA/Monomorphina,Euglenaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Myrmecia,NA/Myrmecia,Trebouxiaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Neocystis,NA/Neocystis,Radiococcaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Nephroselmis,NA,NA/Nephroselmis,Nephroselmidaceae,Nephroselmidales/" plastid.CHECKME.csv
sed -i".bak" "s/Oedogonium,NA/Oedogonium,Oedogoniaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Oltmannsiellopsis,NA/Oltmannsiellopsis,Oltmannsiellopsidaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Pabia,NA/Pabia,Trebouxiaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Parachlorella,NA/Parachlorella,Chlorellaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Paradoxia,NA/Paradoxia,Coccomyxaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Planctonema,NA/Planctonema,Oocystaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Prasinoderma,NA/Prasinoderma,Prasinococcaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Pseudendoclonium,NA/Pseudendoclonium,Kornmanniaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Pseudochloris,NA,NA/Pseudochloris,Chlorellaceae,Chlorellales/" plastid.CHECKME.csv
sed -i".bak" "s/Pyramimonas,NA/Pyramimonas,Pyramimonadaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Stichococcus,NA/Stichococcus,Prasiolaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Stigeoclonium,NA/Stigeoclonium,Chaetophoraceae/" plastid.CHECKME.csv
sed -i".bak" "s/Trachydiscus,NA/Trachydiscus,Pleurochloridaceae/" plastid.CHECKME.csv
sed -i".bak" "s/Watanabea,NA/Watanabea,Trebouxiaceae/" plastid.CHECKME.csv

# Chrysanthemum x morifolium and Pelargonium x hortorum are also weird, but they can be googled.
```

## Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

```text
898 ---------> 894 ---------> 460 ---------> 666
        NA           genus          family
```

```bash
mkdir -p ~/data/organelle/plastid_summary
cd ~/data/organelle/plastid_summary

# filter out accessions without linage information
cat plastid.CHECKME.csv \
    | perl -nl -a -F"," -e \
    '/^#/ and next; ($F[3] eq q{NA} or $F[4] eq q{NA} ) and next; print' \
    > plastid.tmp

# 894
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

# 460
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

# 666
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
# Gossypium herbaceum,2
# Magnolia officinalis,2
# Olea europaea,4
# Oryza sativa,3
# Saccharum hybrid cultivar,2

# strain name not equal to species
cat plastid.DOWNLOAD.csv \
    | grep -v '^#' \
    | perl -nl -a -F"," -e '$F[2] ne $F[3] and print $F[2]' \
    | sort

# Brassica rapa subsp. pekinensis
# Cucumis melo subsp. melo
# Eucalyptus globulus subsp. globulus
# Fagopyrum esculentum subsp. ancestrale
# Fragaria vesca subsp. bracteata
# Fragaria vesca subsp. vesca
# Gossypium herbaceum subsp. africanum
# Hordeum vulgare subsp. vulgare
# Magnolia officinalis subsp. biloba
# Oenothera elata subsp. hookeri
# Olea europaea subsp. cuspidata
# Olea europaea subsp. europaea
# Olea europaea subsp. maroccana
# Olea woodiana subsp. woodiana
# Oryza sativa Indica Group
# Oryza sativa Japonica Group
# Phalaenopsis aphrodite subsp. formosana
# Phyllostachys nigra var. henonis
# Pseudotsuga sinensis var. wilsoniana
# Saccharum hybrid cultivar NCo 310
# Saccharum hybrid cultivar SP80-3280
```

Create abbreviations.

```bash
cd ~/data/organelle/plastid_summary

echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr' > plastid.ABBR.csv
cat plastid.DOWNLOAD.csv \
    | grep -v '^#' \
    | perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 \
    | sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> plastid.ABBR.csv

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

```

## Download sequences and regenerate lineage information.

We don't rename sequences here, so the file has three columns.

And create `plastid_ncbi.csv` with abbr names as taxon file.

```bash
mkdir -p ~/data/organelle/plastid_genomes
cd ~/data/organelle/plastid_genomes

echo "#strain_name,accession,strain_taxon_id" > plastid_name_acc_id.csv
cat ../plastid_summary/plastid.ABBR.csv \
    | grep -v '^#' \
    | perl -nl -a -F"," -e 'print qq{$F[9],$F[1],$F[0]}' \
    | sort \
    >> plastid_name_acc_id.csv

# some warnings fro bioperl, normally just ignore them
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl -f plastid_name_acc_id.csv -p 2>&1 | tee plastid_seq.log

# rsync --progress -av wangq@45.79.80.100:/home/wangq/data/organelle/ ~/data/organelle/

# count downloaded sequences
find . -name "*.fasta" | wc -l

# local, Runtime 10 seconds.
# with --entrez, Runtime 7 minutes and 23 seconds.
# And can't find is still can't find.
cat ../plastid_summary/plastid.ABBR.csv \
    | grep -v '^#' \
    | perl -nl -a -F"," -e 'print qq{$F[0],$F[9]}' \
    | perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file plastid_ncbi.csv
```

## Create alignment plans

We got 47 orders, 58 families, 113 genera, 451 species and **460** accessions.

```bash
cd ~/data/organelle/plastid_summary

# valid genera
cat plastid.ABBR.csv \
    | grep -v "^#" \
    | perl -nl -a -F"," -e \
    '$seen{$F[4]}++; END {for $k (sort keys %seen) {printf qq{,%s,\n}, $k if $seen{$k} > 1}}' \
    > genus.tmp

# intersect between two files
grep -F -f genus.tmp plastid.ABBR.csv > plastid.GENUS.csv

# 460
wc -l plastid.GENUS.csv

#   count every ranks
#   47 order.list.tmp
#   58 family.list.tmp
#  113 genus.list.tmp
#  451 species.list.tmp
cut -d',' -f 4 plastid.GENUS.csv | sort | uniq > species.list.tmp
cut -d',' -f 5 plastid.GENUS.csv | sort | uniq > genus.list.tmp
cut -d',' -f 6 plastid.GENUS.csv | sort | uniq > family.list.tmp
cut -d',' -f 7 plastid.GENUS.csv | sort | uniq > order.list.tmp
wc -l order.list.tmp family.list.tmp genus.list.tmp species.list.tmp

# create again with headers
grep -F -f genus.tmp plastid.ABBR.csv > plastid.GENUS.tmp

# sort by multiply columns, phylum, order, family, abbr
head -n 1 plastid.ABBR.csv > plastid.GENUS.csv
cat plastid.GENUS.tmp \
    | sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> plastid.GENUS.csv

# clean
rm *.tmp *.bak
```

Create alignments without outgroups.

```bash
cd ~/data/organelle/plastid_summary

# tab-seperated
# name  t   qs
cat plastid.GENUS.csv \
    | grep -v "^#" \
    | perl -n -a -F"," -e \
    'BEGIN{ ($g, @s, %h) = (q{}); } chomp for @F; if ($F[4] ne $g) { if ($g) { @s = sort {$h{$a} <=> $h{$b}} @s; $t = shift @s; $qs = join(q{,}, @s); printf qq{%s\t%s\t%s\n}, $g, $t, $qs; } $g = $F[4]; @s = ();} push @s, $F[9]; $h{$F[9]} = $F[0]; END { @s = sort {$h{$a} <=> $h{$b}} @s; $t = shift @s; $qs = join(q{,}, @s); printf qq{%s\t%s\t%s\n}, $g, $t, $qs; }' \
    > genus.tsv

cat plastid.ABBR.csv \
    | grep -v "^#" \
    | perl -n -a -F"," -e \
    'BEGIN{ ($g, @s, %h) = (q{}); } chomp for @F; if ($F[5] ne $g) { if ($g) { @s = sort {$h{$a} <=> $h{$b}} @s; $t = shift @s; $qs = join(q{,}, @s); printf qq{%s\t%s\t%s\n}, $g, $t, $qs; } $g = $F[5]; @s = ();} push @s, $F[9]; $h{$F[9]} = $F[0]; END { @s = sort {$h{$a} <=> $h{$b}} @s; $t = shift @s; $qs = join(q{,}, @s); printf qq{%s\t%s\t%s\n}, $g, $t, $qs; }' \
    > family.tsv

# name  t   qs  o
cat genus.tsv \
    | perl -nl -a -F"\t" -MPath::Tiny -e \
    'BEGIN{ @ls = grep {/\S/} grep {!/^#/} path(q{~/Scripts/withncbi/doc/plastid_OG.md})->lines( { chomp => 1}); for (@ls) {@fs = split(/,/); $h{$fs[0]}= $fs[1];}  } if (exists $h{$F[0]}) { printf qq{%s\t%s\t%s\t%s\n}, $F[0], $F[1], $F[2], $h{$F[0]}; }' \
    > genus_OG.tsv

# every genera
echo -e "mkdir -p ~/data/organelle/plastid.working\ncd ~/data/organelle/plastid.working\n" > ../plastid.cmd.txt
cat genus.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl --seq_dir ~/data/organelle/plastid_genomes --taxon_file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    >> ../plastid.cmd.txt

# this is for finding outgroups
echo -e "mkdir -p ~/data/organelle/plastid_families\ncd ~/data/organelle/plastid_families\n" > ../plastid_families.cmd.txt
cat family.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl --seq_dir ~/data/organelle/plastid_genomes --taxon_file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    >> ../plastid_families.cmd.txt

# genera with outgroups
echo -e "mkdir -p ~/data/organelle/plastid_OG\ncd ~/data/organelle/plastid_OG\n" > ../plastid_OG.cmd.txt
cat genus_OG.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl --seq_dir ~/data/organelle/plastid_genomes --taxon_file ~/data/organelle/plastid_genomes/plastid_ncbi.csv \
    >> ../plastid_OG.cmd.txt

```

## Aligning

### Batch running for genus

The old prepare_run.sh

```bash
mkdir -p ~/data/organelle/plastid.working
cd ~/data/organelle/plastid.working

time sh ../plastid.cmd.txt 2>&1 | tee log_cmd.txt

#----------------------------#
# Approach 1: one by one
#----------------------------#
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo "echo \"====> Processing $d <====\""
    echo sh $d/1_real_chr.sh ; \
    echo sh $d/2_file_rm.sh ; \
    echo sh $d/3_pair_cmd.sh ; \
    echo sh $d/4_rawphylo.sh ; \
    echo sh $d/5_multi_cmd.sh ; \
    echo sh $d/7_multi_db_only.sh ; \
    echo ; \
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

#----------------------------#
# Approach 2: step by step
#----------------------------#
# real_chr
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 1_real_chr.sh | sort `;do \
    echo sh $f ; \
    echo ; \
done  > run_1.sh

# RepeatMasker
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 2_file_rm.sh | sort `;do \
    echo sh $f ; \
    echo ; \
done  > run_2.sh

# pair
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 3_pair_cmd.sh | sort `;do \
    echo sh $f ; \
    echo ; \
done  > run_3.sh

# rawphylo
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 4_rawphylo.sh | sort `;do \
    echo sh $f ; \
    echo ; \
done  > run_4.sh

# multi cmd
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 5_multi_cmd.sh | sort `;do \
    echo sh $f ; \
    echo ; \
done  > run_5.sh

# multi db
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 7_multi_db_only.sh | sort `;do \
    echo sh $f ; \
    echo ; \
done  > run_7.sh

cat run_1.sh | grep . | parallel -j 4 2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel -j 2 2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel -j 1 2>&1 | tee log_3.txt
cat run_4.sh | grep . | parallel -j 1 2>&1 | tee log_4.txt
cat run_5.sh | grep . | parallel -j 1 2>&1 | tee log_5.txt
cat run_7.sh | grep . | parallel -j 2 2>&1 | tee log_7.txt

#----------------------------#
# Charting on Windows
#----------------------------#
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do \
    export d_base=`basename $d` ; \
    echo "perl d:/Scripts/fig_table/collect_common_basic.pl    -d $d_base" ; \
    echo "perl d:/Scripts/alignDB/stat/common_chart_factory.pl -i $d_base/$d_base.common.xlsx" ; \
    echo "perl d:/Scripts/alignDB/stat/multi_chart_factory.pl  -i $d_base/$d_base.multi.xlsx" ; \
    echo "perl d:/Scripts/alignDB/stat/gc_chart_factory.pl     -i $d_base/$d_base.gc.xlsx" ; \
    echo ; \
done  > run_chart.bat
perl -pi -e 's/\n/\r\n/g' run_chart.bat

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | xargs rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | xargs rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | xargs rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | xargs rm
```

### Self alignments.

```bash
cd ~/data/organelle/

perl -p -e 's/plastid\.working/plastid_self.working/g; s/strain_bz/strain_bz_self/g; s/(\-\-use_name)/\1 --length 1000 /g;' plastid.cmd.txt > plastid_self.cmd.txt

mkdir -p ~/data/organelle/plastid_self.working
cd ~/data/organelle/plastid_self.working

time sh ../plastid_self.cmd.txt 2>&1 | tee log_cmd.txt

# Don't need 6_feature_cmd.sh
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo "echo \"====> Processing $d <====\""
    echo sh $d/1_real_chr.sh ; \
    echo sh $d/2_file_rm.sh ; \
    echo sh $d/3_self_cmd.sh ; \
    echo sh $d/4_proc_cmd.sh ; \
    echo sh $d/5_circos_cmd ; \
    echo sh $d/7_pair_stat.sh ; \
    echo ; \
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

#----------------------------#
# Charting on Windows
#----------------------------#
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do \
    export d_base=`basename $d` ; \
    echo "perl d:/Scripts/fig_table/collect_common_basic.pl    -d $d_base" ; \
    echo "perl d:/Scripts/alignDB/stat/common_chart_factory.pl -i $d_base/${d_base}_paralog.common.xlsx" ; \
    echo "perl d:/Scripts/alignDB/stat/gc_chart_factory.pl     -i $d_base/${d_base}_paralog.gc.xlsx" ; \
    echo ; \
done  > run_chart.bat
perl -pi -e 's/\n/\r\n/g' run_chart.bat

# clean
find . -mindepth 1 -maxdepth 2 -type d -name "*_raw" | xargs rm -fr
find . -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | xargs rm -fr

```

### Alignments of families for outgroups.

The old prepare_run.sh

```bash
mkdir -p ~/data/organelle/plastid_families
cd ~/data/organelle/plastid_families

time sh ../plastid_families.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo "echo \"====> Processing $d <====\""
    echo sh $d/1_real_chr.sh ; \
    echo sh $d/2_file_rm.sh ; \
    echo sh $d/3_pair_cmd.sh ; \
    echo sh $d/4_rawphylo.sh ; \
    echo sh $d/5_multi_cmd.sh ; \
    echo ; \
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do \
    export d_base=`basename $d` ; \
    echo "perl d:/Scripts/fig_table/collect_common_basic.pl    -d $d_base" ; \
    echo ; \
done  > run_chart.bat
perl -pi -e 's/\n/\r\n/g' run_chart.bat

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    export d_base=`basename $d` ;
    export f_base="${d_base}/${d_base}_phylo/${d_base}" ;
    if [ -f $f_base.nwk ]
    then
        echo $f_base ;  
        nw_display -s -b 'visibility:hidden' $f_base.nwk > $f_base.svg ;
    fi
done

find . -type f -path "*_phylo*" -name "*.nwk"

```

After manually editing.

*D* of outgroups are around 0.05.

```bash
mkdir -p ~/data/organelle/plastid_OG
cd ~/data/organelle/plastid_OG

time sh ../plastid_OG.cmd.sh 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `; do
    echo "echo \"====> Processing $d <====\""
    echo sh $d/1_real_chr.sh ;
    echo sh $d/2_file_rm.sh ;
    echo sh $d/3_pair_cmd.sh ;
    echo sh $d/4_rawphylo.sh ;
    echo sh $d/5_multi_cmd.sh ;
    echo sh $d/7_multi_db_only.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

```

### summary

Create `plastid.list.csv` from `plastid.GENUS.csv` with sequence lengths.

```bash
mkdir -p ~/data/organelle/plastid_summary
cd ~/data/organelle/plastid_summary

find ~/data/organelle/plastid.working -type f -name "chr.sizes" | sort \
    | xargs perl -nl -e 'BEGIN{print q{genus,strain_abbr,accession,length}}; $_ =~ s/\t/\,/; $ARGV =~ /working\/(\w+)\/(\w+)\//; print qq{$1,$2,$_}' > length.tmp

perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t ~/data/organelle/plastid_genomes/plastid.GENUS.csv -m length.tmp -f 1 -f2 2 --concat --stdout \
    | perl -nl -a -F"," -e 'print qq{$F[5],$F[4],$F[2],$F[0],$F[1],$F[13]}' \
    >  plastid.list.csv

```
