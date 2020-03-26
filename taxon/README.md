# Process plastid genomes

[TOC levels=1-3]: # ""

- [Process plastid genomes](#process-plastid-genomes)
- [Work flow.](#work-flow)
- [Update taxdmp](#update-taxdmp)
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
- [Add lineage information](#add-lineage-information)
  - [Can't get clear taxon information](#cant-get-clear-taxon-information)
- [Filtering based on valid families and genera](#filtering-based-on-valid-families-and-genera)
- [Find a way to name these](#find-a-way-to-name-these)
- [Download sequences and regenerate lineage information](#download-sequences-and-regenerate-lineage-information)
- [Prepare sequences for lastz](#prepare-sequences-for-lastz)
- [Aligning without outgroups](#aligning-without-outgroups)
  - [Create alignments plans without outgroups](#create-alignments-plans-without-outgroups)
  - [Batch running for genera](#batch-running-for-genera)
  - [Alignments of families for outgroups.](#alignments-of-families-for-outgroups)
- [Aligning with outgroups](#aligning-with-outgroups)
  - [Create `plastid_OG.md` for picking outgroups](#create-plastid_ogmd-for-picking-outgroups)
  - [Create alignments plans with outgroups](#create-alignments-plans-with-outgroups)
- [Self alignments](#self-alignments)
- [LSC and SSC](#lsc-and-ssc)
  - [Can't get clear IR information](#cant-get-clear-ir-information)
  - [Slices of IR, LSC and SSC](#slices-of-ir-lsc-and-ssc)
- [Cyanobacteria](#cyanobacteria)
  - [Genus and Species counts](#genus-and-species-counts)
- [Summary](#summary)
  - [Copy xlsx files](#copy-xlsx-files)
  - [Genome list](#genome-list)
  - [Statistics of genome alignments](#statistics-of-genome-alignments)
  - [Groups](#groups)
  - [Phylogenic trees of each genus with outgroup](#phylogenic-trees-of-each-genus-with-outgroup)
  - [d1, d2](#d1-d2)


The following command lines are about how I processed the plastid genomes of green plants. Many
tools of `taxon/` are used here, which makes a good example for users.

# Work flow.

```text
id ---> lineage ---> filtering ---> naming ---> strain_info.pl   ---> egaz/multi_batch.pl
                                      |                                 ^
                                      |-------> batch_get_seq.pl -------|
```

I'm sure there are no commas in names. So for convenient, don't use Text::CSV_XS.

# Update taxdmp

*Update `~/data/NCBI/taxdmp` before running `id_restrict.pl` or `id_project_to.pl`*.

# Scrap id and acc from NCBI

Open browser and visit
[NCBI plastid page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=plastid).
Save page to a local file, html only. In this case, it's `doc/green_plants_plastid_200326.html`.

All [Eukaryota](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?opt=plastid&taxid=2759),
`doc/eukaryota_plastid_200326.html`.

Or [this link](http://www.ncbi.nlm.nih.gov/genome/browse/?report=5).

```text
Eukaryota (2759)                4711
    Viridiplantae (33090)       4437
        Chlorophyta (3041)      151
        Streptophyta (35493)    4284
```

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions from all history
pages.

```csv
id,acc
996148,NC_017006
```

Got **4755** accessions.

```bash
mkdir -p ~/data/plastid/GENOMES
cd ~/data/plastid/GENOMES

rm webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/eukaryota_plastid_200326.html \
    >> webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/green_plants_plastid_200326.html \
    >> webpage_id_seq.csv

```

Use `taxon/gb_taxon_locus.pl` to extract information from refseq genbank files.

```bash
cd ~/data/plastid/GENOMES

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.1.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.2.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.3.genomic.gbff.gz

gzip -dcf plastid.*.genomic.gbff.gz > genomic.gbff

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv

rm genomic.gbff
#rm *.genomic.gbff.gz

cat refseq_id_seq.csv | grep -v "^#" | wc -l
# 4722

# combine
cat webpage_id_seq.csv refseq_id_seq.csv |
    sort -u | # duplicated id-seq pair
    sort -t, -k1,1 \
    > id_seq.csv

cat id_seq.csv | grep -v "^#" | wc -l
# 4755

```

# Add lineage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `CHECKME.csv`.

```bash
mkdir -p ~/data/plastid/summary
cd ~/data/plastid/summary

# generate a .csv file for manually checking
echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum' > CHECKME.csv
cat ../GENOMES/id_seq.csv |
    grep -v "^#" |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank species |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank genus |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank family |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank order |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank class |
    perl ~/Scripts/withncbi/taxon/id_project_to.pl -s "," --rank phylum |
    sort -t',' -k9,9 -k8,8 -k7,7 -k6,6 -k5,5 \
    >> CHECKME.csv

```

Manually correct lineages.

Taxonomy information from [AlgaeBase](http://www.algaebase.org),
[Wikipedia](https://www.wikipedia.org/) and [Encyclopedia of Life](http://eol.org/).

```bash
cd ~/data/plastid/summary

# darwin (bsd) need "" for -i
sed -i".bak" "s/\'//g" CHECKME.csv

# Koliella corcontica (a green algae) was grouped to Streptophyta.
# Koliella longiseta
perl -pi -e 's/Koliella,\w+,\w+,\w+,\w+/Koliella,Klebsormidiaceae,Klebsormidiales,Klebsormidiophyceae,Chlorophyta/g' CHECKME.csv
sed -i".bak" "s/Klebsormidiophyceae,Streptophyta/Klebsormidiophyceae,Chlorophyta/" CHECKME.csv

# Chrysanthemum x morifolium and Pelargonium x hortorum are also weird, but they can be googled.

# missing families
# queried from http://www.algaebase.org/
sed -i".bak" "s/Aureococcus,NA/Aureococcus,Pelagomonadaceae/" CHECKME.csv
sed -i".bak" "s/Aureoumbra,NA/Aureoumbra,Sarcinochrysidaceae/" CHECKME.csv
sed -i".bak" "s/Bigelowiella,NA/Bigelowiella,Chlorarachniaceae/" CHECKME.csv
sed -i".bak" "s/Choricystis,NA/Choricystis,Coccomyxaceae/" CHECKME.csv
sed -i".bak" "s/Carteria,NA,NA/Carteria,Chlamydomonadaceae,Chlamydomonadales/" CHECKME.csv
sed -i".bak" "s/Cryptoglena,NA,NA/Cryptoglena,Euglenaceae,Euglenales/" CHECKME.csv
sed -i".bak" "s/Cymbomonas,NA/Cymbomonas,Pyramimonadaceae/" CHECKME.csv
sed -i".bak" "s/Dicloster,NA,NA/Dicloster,Chlorellaceae,Chlorellales/" CHECKME.csv
sed -i".bak" "s/Dictyochloropsis,NA/Dictyochloropsis,Trebouxiaceae/" CHECKME.csv
sed -i".bak" "s/Euglenaformis,NA/Euglenaformis,Euglenaceae/" CHECKME.csv
sed -i".bak" "s/Euglenaria,NA/Euglenaria,Euglenaceae/" CHECKME.csv
sed -i".bak" "s/Eutreptiella,NA/Eutreptiella,Eutreptiaceae/" CHECKME.csv
sed -i".bak" "s/Fusochloris,NA/Fusochloris,Microthamniaceae/" CHECKME.csv
sed -i".bak" "s/Geminella,NA,NA/Geminella,Chlorellaceae,Chlorellales/" CHECKME.csv
sed -i".bak" "s/Gloeotilopsis,NA/Gloeotilopsis,Ulotrichaceae/" CHECKME.csv
sed -i".bak" "s/Helicosporidium,NA/Helicosporidium,Chlorellaceae/" CHECKME.csv
sed -i".bak" "s/Microthamnion,NA/Microthamnion,Microthamniaceae/" CHECKME.csv
sed -i".bak" "s/Monomorphina,NA/Monomorphina,Euglenaceae/" CHECKME.csv
sed -i".bak" "s/Myrmecia,NA/Myrmecia,Trebouxiaceae/" CHECKME.csv
sed -i".bak" "s/Neocystis,NA/Neocystis,Radiococcaceae/" CHECKME.csv
sed -i".bak" "s/Nephroselmis,NA,NA/Nephroselmis,Nephroselmidaceae,Nephroselmidales/" CHECKME.csv
sed -i".bak" "s/Oedogonium,NA/Oedogonium,Oedogoniaceae/" CHECKME.csv
sed -i".bak" "s/Oedocladium,NA/Oedocladium,Oedogoniaceae/" CHECKME.csv
sed -i".bak" "s/Oltmannsiellopsis,NA/Oltmannsiellopsis,Oltmannsiellopsidaceae/" CHECKME.csv
sed -i".bak" "s/Pabia,NA/Pabia,Trebouxiaceae/" CHECKME.csv
sed -i".bak" "s/Parachlorella,NA/Parachlorella,Chlorellaceae/" CHECKME.csv
sed -i".bak" "s/Paradoxia,NA/Paradoxia,Coccomyxaceae/" CHECKME.csv
sed -i".bak" "s/Planctonema,NA/Planctonema,Oocystaceae/" CHECKME.csv
sed -i".bak" "s/Prasinoderma,NA/Prasinoderma,Prasinococcaceae/" CHECKME.csv
sed -i".bak" "s/Pseudendoclonium,NA/Pseudendoclonium,Kornmanniaceae/" CHECKME.csv
sed -i".bak" "s/Pseudochloris,NA,NA/Pseudochloris,Chlorellaceae,Chlorellales/" CHECKME.csv
sed -i".bak" "s/Pyramimonas,NA/Pyramimonas,Pyramimonadaceae/" CHECKME.csv
sed -i".bak" "s/Stichococcus,NA/Stichococcus,Prasiolaceae/" CHECKME.csv
sed -i".bak" "s/Stigeoclonium,NA/Stigeoclonium,Chaetophoraceae/" CHECKME.csv
sed -i".bak" "s/Trachydiscus,NA/Trachydiscus,Pleurochloridaceae/" CHECKME.csv
sed -i".bak" "s/Verdigellas,NA/Verdigellas,Palmophyllaceae/" CHECKME.csv
sed -i".bak" "s/Watanabea,NA/Watanabea,Trebouxiaceae/" CHECKME.csv

sed -i".bak" "s/Chlorellales,NA/Chlorellales,Trebouxiophyceae/" CHECKME.csv
sed -i".bak" "s/Ettlia,NA/Ettlia,Chlorococcaceae/" CHECKME.csv

# Chlorophyceae, incertae sedis
#sed -i".bak" "s/Jenufa,NA,NA/Jenufa,NA,NA/" CHECKME.csv
# Chlorophyta incertae sedis
#sed -i".bak" "s/Picocystis,NA,NA/Picocystis,NA,NA/" CHECKME.csv
sed -i".bak" "s/Pleurastrum,NA,NA/Pleurastrum,NA,Chlamydomonadales/" CHECKME.csv

# missing orders
sed -i".bak" "s/Leptocylindraceae,NA/Leptocylindraceae,Chaetocerotales/" CHECKME.csv
sed -i".bak" "s/Rhizosoleniaceae,NA/Rhizosoleniaceae,Rhizosoleniales/" CHECKME.csv
sed -i".bak" "s/Babesiidae,NA/Babesiidae,Piroplasmida/" CHECKME.csv
sed -i".bak" "s/Theileriidae,NA/Theileriidae,Piroplasmida/" CHECKME.csv
sed -i".bak" "s/Treubariaceae,NA/Treubariaceae,Chlorococcales/" CHECKME.csv
sed -i".bak" "s/Oltmannsiellopsidaceae,NA/Oltmannsiellopsidaceae,Oltmannsiellopsidales/" CHECKME.csv
sed -i".bak" "s/Pycnococcaceae,NA/Pycnococcaceae,Pseudoscourfieldiales/" CHECKME.csv
sed -i".bak" "s/Coccomyxaceae,NA/Coccomyxaceae,Chlorococcales/" CHECKME.csv

# missing classes and phylums
sed -i".bak" "s/Glaucocystophyceae,NA/Glaucocystophyceae,Glaucophyta/" CHECKME.csv
sed -i".bak" "s/Dinophyceae,NA/Dinophyceae,Dinoflagellata/" CHECKME.csv

sed -i".bak" "s/Apiales,NA/Apiales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Aquifoliales,NA/Aquifoliales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Asterales,NA/Asterales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Austrobaileyales,NA/Austrobaileyales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Brassicales,NA/Brassicales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Buxales,NA/Buxales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Canellales,NA/Canellales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Caryophyllales,NA/Caryophyllales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Celastrales,NA/Celastrales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Cornales,NA/Cornales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Cucurbitales,NA/Cucurbitales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Dipsacales,NA/Dipsacales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Ericales,NA/Ericales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Fabales,NA/Fabales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Gentianales,NA/Gentianales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Geraniales,NA/Geraniales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Ginkgoales,NA/Ginkgoales,Ginkgoopsida/" CHECKME.csv
sed -i".bak" "s/Gnetales,NA/Gnetales,Gnetopsida/" CHECKME.csv
sed -i".bak" "s/Gleicheniales,NA/Gleicheniales,Gleicheniales/" CHECKME.csv
sed -i".bak" "s/Isoetales,NA/Isoetales,Isoetopsida/" CHECKME.csv
sed -i".bak" "s/Lamiales,NA/Lamiales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Lycopodiales,NA/Lycopodiales,Lycopodiopsida/" CHECKME.csv
sed -i".bak" "s/Magnoliales,NA/Magnoliales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Malpighiales,NA/Malpighiales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Malvales,NA/Malvales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Marattiales,NA/Marattiales,Marattiopsida/" CHECKME.csv
sed -i".bak" "s/Myrtales,NA/Myrtales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Ophioglossales,NA/Ophioglossales,Ophioglossales/" CHECKME.csv
sed -i".bak" "s/Osmundales,NA/Osmundales,Polypodiopsida/" CHECKME.csv
sed -i".bak" "s/Pinales,NA/Pinales,Pinopsida/" CHECKME.csv
sed -i".bak" "s/Polypodiales,NA/Polypodiales,Polypodiopsida/" CHECKME.csv
sed -i".bak" "s/Proteales,NA/Proteales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Ranunculales,NA/Ranunculales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Rosales,NA/Rosales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Salviniales,NA/Salviniales,Polypodiopsida/" CHECKME.csv
sed -i".bak" "s/Saxifragales,NA/Saxifragales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Schizaeales,NA/Schizaeales,Polypodiopsida/" CHECKME.csv
sed -i".bak" "s/Selaginellales,NA/Selaginellales,Isoetopsida/" CHECKME.csv
sed -i".bak" "s/Solanales,NA/Solanales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Trochodendrales,NA/Trochodendrales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Vitales,NA/Vitales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Welwitschiales,NA/Welwitschiales,Gnetopsida/" CHECKME.csv
sed -i".bak" "s/Zygophyllales,NA/Zygophyllales,Magnoliopsida/" CHECKME.csv
sed -i".bak" "s/Araucariales,NA/Araucariales,Pinidae/g" CHECKME.csv
sed -i".bak" "s/Cupressales,NA/Cupressales,Pinidae/g" CHECKME.csv
sed -i".bak" "s/Cycadales,NA/Cycadales,Cycadopsida/g" CHECKME.csv
sed -i".bak" "s/Ephedrales,NA/Ephedrales,Gnetidae/g" CHECKME.csv
sed -i".bak" "s/Ceratophyllales,NA/Ceratophyllales,Magnoliopsida/g" CHECKME.csv
sed -i".bak" "s/Chloranthales.NA/Chloranthales,Magnoliopsida/g" CHECKME.csv
sed -i".bak" "s/Fagales,NA/Fagales,Eudicotyledoneae/g" CHECKME.csv
sed -i".bak" "s/Garryales,NA/Garryales,Magnoliopsida/g" CHECKME.csv
sed -i".bak" "s/Huerteales,NA/Huerteales,Magnoliopsida/g" CHECKME.csv
sed -i".bak" "s/Icacinales,NA/Icacinales,eudicots/g" CHECKME.csv
sed -i".bak" "s/Laurales,NA/Laurales,Magnoliopsida/g" CHECKME.csv
sed -i".bak" "s/Nymphaeales,NA/Nymphaeales,Magnoliopsida/g" CHECKME.csv
sed -i".bak" "s/Oxalidales,NA/Oxalidales,eudicots/g" CHECKME.csv
sed -i".bak" "s/Piperales,NA/Piperales,Magnoliopsida/g" CHECKME.csv
sed -i".bak" "s/Santalales,NA/Santalales,eudicots/g" CHECKME.csv
sed -i".bak" "s/Sapindales,NA/Sapindales,eudicots/g" CHECKME.csv
sed -i".bak" "s/Amborellales,NA/Amborellales,Magnoliopsida/g" CHECKME.csv

# Entry Merged. Taxid 1605147 was merged into taxid 142389 on October 16, 2015.
sed -i".bak" "/1605147,/d" CHECKME.csv

# comma in names
sed -i".bak" "/167339,/d" CHECKME.csv

# missing all
sed -i".bak" "/2003521,/d" CHECKME.csv

# Cercozoa 丝足虫门
sed -i".bak" "s/Chlorarachniophyceae,NA/Chlorarachniophyceae,Cercozoa/" CHECKME.csv

sed -i".bak" "s/Euglyphida,NA,NA/Euglyphida,Filosa,Cercozoa/" CHECKME.csv
sed -i".bak" "s/Gymnochlora,NA,NA,NA,NA/Gymnochlora,Chlorarachniaceae,Chlorarachniales,Chlorarachniophyceae,Cercozoa/" CHECKME.csv
sed -i".bak" "s/Lotharella,NA,NA,NA,NA/Lotharella,Chlorarachniaceae,Chlorarachniales,Chlorarachniophyceae,Cercozoa/" CHECKME.csv
sed -i".bak" "s/Partenskyella,NA,NA,NA,NA/Partenskyella,Chlorarachniaceae,Chlorarachniales,Chlorarachniophyceae,Cercozoa/" CHECKME.csv
sed -i".bak" "s/Bigelowiella,Chlorarachniaceae,NA,NA,NA/Bigelowiella,Chlorarachniaceae,Chlorarachniales,Chlorarachniophyceae,Cercozoa/" CHECKME.csv

# Haptophyta 定鞭藻门
sed -i".bak" "s/Haptophyceae,NA/Haptophyceae,Haptophyta/" CHECKME.csv
sed -i".bak" "s/Isochrysidales,NA,NA/Isochrysidales,Coccolithophyceae,Haptophyta/" CHECKME.csv
sed -i".bak" "s/Pavlovales,NA,NA/Pavlovales,Pavlovophyceae,Haptophyta/" CHECKME.csv
sed -i".bak" "s/Phaeocystales,NA,NA/Phaeocystales,Coccolithophyceae,Haptophyta/" CHECKME.csv

# Ochrophyta 褐藻门
parallel -j 1 '
    sed -i".bak" "s/{},NA/{},Ochrophyta/" CHECKME.csv
    ' ::: \
        Bolidophyceae Dictyochophyceae Eustigmatophyceae \
        Pelagophyceae Phaeophyceae Raphidophyceae \
        Synurophyceae Xanthophyceae

# Rhodophyta 红藻门
parallel -j 1 '
    sed -i".bak" "s/{},NA/{},Rhodophyta/" CHECKME.csv
    ' ::: \
        Bangiophyceae Compsopogonophyceae Florideophyceae \
        Rhodellophyceae Stylonematophyceae

# Cryptophyta 隐藻门
parallel -j 1 '
    sed -i".bak" "s/{},NA/{},Cryptophyta/" CHECKME.csv
    ' ::: \
        Cryptophyta Cryptophyceae
sed -i".bak" "s/Cryptomonadales,NA,NA/Cryptomonadales,Cryptophyceae,Cryptophyta/" CHECKME.csv
sed -i".bak" "s/Pyrenomonadales,NA,NA/Pyrenomonadales,Cryptophyceae,Cryptophyta/" CHECKME.csv

# Charophyta 轮藻门
sed -i".bak" "s/Charophyceae,Streptophyta/Charophyceae,Charophyta/" CHECKME.csv
sed -i".bak" "s/Chlorokybophyceae,Streptophyta/Chlorokybophyceae,Charophyta/" CHECKME.csv
sed -i".bak" "s/Coleochaetophyceae,Streptophyta/Coleochaetophyceae,Charophyta/" CHECKME.csv
sed -i".bak" "s/Zygnemophyceae,Streptophyta/Zygnemophyceae,Charophyta/" CHECKME.csv

# Chlorophyta 绿藻门
sed -i".bak" "s/Mesostigmatophyceae,Streptophyta/Mesostigmatophyceae,Chlorophyta/" CHECKME.csv

```

Split Streptophyta according to classical plant classification.

* Streptophyta
  * Streptophytina
    * Embryophyta
      * Tracheophyta
        * Euphyllophyta
          * Spermatophyta
            * Magnoliopsida (flowering plants) 3398 - Angiosperm
            * Acrogymnospermae 1437180 - Gymnosperm
          * Polypodiopsida 241806 - ferns
        * Lycopodiopsida 1521260 - clubmosses
      * Anthocerotophyta 13809 - hornworts
      * Bryophyta 3208 - mosses
      * Marchantiophyta 3195 - liverworts

```bash
cd ~/data/plastid/summary

# Angiosperms
perl ~/Scripts/withncbi/taxon/id_members.pl 3398 --rank family |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Angiosperms/g
        '\'' CHECKME.csv
    '

# Gymnosperms
perl ~/Scripts/withncbi/taxon/id_members.pl 1437180 --rank family |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Gymnosperms/g
        '\'' CHECKME.csv
    '

# Pteridophytes
perl ~/Scripts/withncbi/taxon/id_members.pl 241806 --rank family |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Pteridophytes/g
        '\'' CHECKME.csv
    '

perl ~/Scripts/withncbi/taxon/id_members.pl 1521260 --rank family |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Pteridophytes/g
        '\'' CHECKME.csv
    '

# Bryophytes
perl ~/Scripts/withncbi/taxon/id_members.pl 13809 --rank family |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Bryophytes/g
        '\'' CHECKME.csv
    '

perl ~/Scripts/withncbi/taxon/id_members.pl 3208 --rank family |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Bryophytes/g
        '\'' CHECKME.csv
    '

perl ~/Scripts/withncbi/taxon/id_members.pl 3195 --rank family |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Bryophytes/g
        '\'' CHECKME.csv
    '

rm *.tmp *.bak

```

## Can't get clear taxon information

* Genus
  + Chromera
  + Elliptochloris
  + Ettlia
  + Picocystis
  + Xylochloris
  + Jenufa
  + Pleurastrum

* Species
  + Chromerida sp. RM11
  + Trebouxiophyceae sp. MX-AZ01
  + Trebouxiophyceae sp. TP-2016a

# Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

```text
4754 ---------> 4733 ---------> 3338 ---------> 4292
        NA             genus          family
```

```bash
mkdir -p ~/data/plastid/summary
cd ~/data/plastid/summary

cat CHECKME.csv | grep -v "^#" | wc -l
# 4754

# filter out accessions without linage information (strain, species, genus and family)
cat CHECKME.csv |
    perl -nla -F"," -e '
        /^#/ and next;
        ($F[2] eq q{NA} or $F[3] eq q{NA} or $F[4] eq q{NA} or $F[5] eq q{NA} ) and next;
        print
    ' \
    > valid.tmp

wc -l valid.tmp
# 4733

#----------------------------#
# Genus
#----------------------------#
# valid genera
cat valid.tmp |
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
grep -F -f genus.tmp valid.tmp > valid.genus.tmp

wc -l valid.genus.tmp
# 3338

#----------------------------#
# Family
#----------------------------#
# get some genera back as candidates for outgroup
cat valid.genus.tmp |
    perl -nla -F"," -e 'printf qq{,$F[5],\n}' \
    > family.tmp

# intersect between two files
grep -F -f family.tmp valid.tmp > valid.family.tmp

wc -l valid.family.tmp
# 4292

#----------------------------#
# results produced in this step
#----------------------------#
head -n 1 CHECKME.csv > DOWNLOAD.csv
cat valid.family.tmp >> DOWNLOAD.csv

# clean
rm *.tmp *.bak

```

# Find a way to name these

Seems it's OK to use species as names.

```bash
cd ~/data/plastid/summary

# sub-species
cat DOWNLOAD.csv |
    perl -nl -a -F"," -e '
        /^#/i and next;
        $seen{$F[3]}++;
        END {
            for $k (keys %seen){printf qq{%s,%d\n}, $k, $seen{$k} if $seen{$k} > 1}
        };
    ' |
    sort
#Allium cyathophorum,2
#Anemone hepatica,2
#Arabidopsis lyrata,2
#Astragalus mongholicus,2
#Brassica rapa,2
#Cannabis sativa,2
#Capsicum baccatum,3
#Changiostyrax dolichocarpus,2
#Corallorhiza striata,2
#Corylus ferox,2
#Fragaria vesca,2
#Gossypium herbaceum,2
#Hippophae rhamnoides,2
#Hordeum vulgare,2
#Magnolia officinalis,2
#Marchantia polymorpha,2
#Musa balbisiana,2
#Olea europaea,4
#Oryza sativa,4
#Panax japonicus,2
#Paris polyphylla,2
#Physcomitrella patens,2
#Pisum sativum,2
#Plasmodium falciparum,2
#Saccharum hybrid cultivar,3
#Sanguisorba tenuifolia,2
#Sinalliaria limprichtiana,2
#Solanum bukasovii,2
#Solanum lycopersicum,2
#Solanum stenotomum,2
#Sophora alopecuroides,2
#Vitis aestivalis,2
#Vitis cinerea,4
#Vitis rotundifolia,2
#Wurfbainia villosa,2

# strain name not equal to species
cat DOWNLOAD.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e '$F[2] ne $F[3] and print $F[2]' |
    sort
#Actinidia callosa var. henryi
#Allium cyathophorum var. farreri
#Anemone cernua var. koreana
#Anemone hepatica var. asiatica
#Anemone hepatica var. japonica
#Arabidopsis lyrata subsp. lyrata
#Astragalus mongholicus var. nakaianus
#Babesia bovis T2Bo
#Babesia microti strain RI
#Betula pendula var. carelica
#Brassica rapa subsp. pekinensis
#Calycanthus floridus var. glaucus
#Calypso bulbosa var. occidentalis
#Capsicum baccatum var. baccatum
#Capsicum baccatum var. pendulum
#Capsicum baccatum var. praetermissum
#Caragana rosea var. rosea
#Corallorhiza striata var. involuta
#Corallorhiza striata var. striata
#Corylus ferox var. thibetica
#Cucumis melo subsp. melo
#Eucalyptus globulus subsp. globulus
#Fagopyrum esculentum subsp. ancestrale
#Fragaria vesca subsp. bracteata
#Fragaria vesca subsp. vesca
#Gossypium herbaceum subsp. africanum
#Gracilaria tenuistipitata var. liui
#Hippophae rhamnoides subsp. yunnanensis
#Hordeum vulgare subsp. spontaneum
#Hordeum vulgare subsp. vulgare
#Lilium martagon var. pilosiusculum
#Magnolia macrophylla var. dealbata
#Magnolia officinalis subsp. biloba
#Marchantia polymorpha subsp. ruderalis
#Musa balbisiana var. balbisiana
#Oenothera elata subsp. hookeri
#Olea europaea subsp. cuspidata
#Olea europaea subsp. europaea
#Olea europaea subsp. maroccana
#Olea woodiana subsp. woodiana
#Oryza sativa f. spontanea
#Oryza sativa Indica Group
#Oryza sativa Japonica Group
#Panax japonicus var. bipinnatifidus
#Paris polyphylla var. chinensis
#Paris polyphylla var. yunnanensis
#Phalaenopsis aphrodite subsp. formosana
#Phryma leptostachya subsp. asiatica
#Phyllostachys nigra var. henonis
#Pisum sativum subsp. elatius
#Plasmodium chabaudi chabaudi
#Plasmodium falciparum 3D7
#Plasmodium falciparum HB3
#Pseudotsuga sinensis var. wilsoniana
#Rosa chinensis var. spontanea
#Saccharum hybrid cultivar NCo 310
#Saccharum hybrid cultivar SP80-3280
#Sanguisorba tenuifolia var. alba
#Sinalliaria limprichtiana var. grandifolia
#Solanum bukasovii f. multidissectum
#Solanum stenotomum subsp. goniocalyx
#Sophora alopecuroides var. alopecuroides
#Styphnolobium japonicum var. japonicum
#Thalassiosira oceanica CCMP1005
#Trichopus zeylanicus subsp. travancoricus
#Vitis aestivalis var. linsecomii
#Vitis cinerea var. cinerea
#Vitis cinerea var. floridana
#Vitis cinerea var. helleri
#Vitis rotundifolia var. munsoniana
#Wurfbainia villosa var. xanthioides

```

Create abbreviations.

```bash
cd ~/data/plastid/summary

echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr' > ABBR.csv
cat DOWNLOAD.csv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 --shortsub |
    sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> ABBR.csv

```

# Download sequences and regenerate lineage information

We don't rename sequences here, so the file has three columns.

And create `taxon_ncbi.csv` with abbr names as taxon file.

```bash
cd ~/data/organelle/plastid/GENOMES

echo "#strain_name,accession,strain_taxon_id" > name_acc_id.csv
cat ../summary/ABBR.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e 'print qq{$F[9],$F[1],$F[0]}' |
    sort \
    >> name_acc_id.csv

# local, Runtime 10 seconds.
# with --entrez, Runtime 7 minutes and 23 seconds.
# And which-can't-find is still which-can't-find.
cat ../summary/ABBR.csv |
    grep -v '^#' |
    perl -nla -F"," -e 'print qq{$F[0],$F[9]}' |
    uniq |
    perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file taxon_ncbi.csv

# Some warnings about trans-splicing genes from BioPerl, just ignore them
# eutils restricts 3 connections
cat name_acc_id.csv |
    grep -v '^#' |
    2>&1 parallel --colsep ',' --no-run-if-empty --linebuffer -k -j 3 "
        echo -e '==> id: [{1}]\tseq: [{2}]\n'
        mkdir -p {1}
        if [[ -e '{1}/{2}.gff' && -e '{1}/{2}.fa' ]] ; then
            echo -e '    Sequence [{1}/{2}] exists, next\n'
            exit
        fi
        
        # gb
        echo -e '    [{1}/{2}].gb'
        curl -Ls \
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={2}&rettype=gb&retmode=text' \
            > {1}/{2}.gb
        
        # fasta
        echo -e '    [{1}/{2}].fa'
        curl -Ls \
            'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={2}&rettype=fasta&retmode=text' \
            > {1}/{2}.fa
        
        # gff
        echo -e '    [{1}/{2}].gff'
        perl ~/Scripts/withncbi/taxon/bp_genbank2gff3.pl {1}/{2}.gb -o stdout > {1}/{2}.gff
        perl -i -nlp -e '/^\#\#FASTA/ and last' {1}/{2}.gff
        
        echo
    " |
    tee download_seq.log

# count downloaded sequences
find . -name "*.fa" | wc -l

```

Numbers for higher ranks are: 84 orders, 165 families, 466 genera and 2225 species.

```bash
cd ~/data/organelle/plastid/summary

# valid genera
cat ABBR.csv |
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
grep -F -f genus.tmp ABBR.csv > GENUS.csv

# 2253
wc -l GENUS.csv

# count every ranks
#   84 order.list.tmp
#  165 family.list.tmp
#  466 genus.list.tmp
# 2225 species.list.tmp
cut -d',' -f 4 GENUS.csv | sort | uniq > species.list.tmp
cut -d',' -f 5 GENUS.csv | sort | uniq > genus.list.tmp
cut -d',' -f 6 GENUS.csv | sort | uniq > family.list.tmp
cut -d',' -f 7 GENUS.csv | sort | uniq > order.list.tmp
wc -l order.list.tmp family.list.tmp genus.list.tmp species.list.tmp

# create again with headers
grep -F -f genus.tmp ABBR.csv > GENUS.tmp

# sort by multiply columns, phylum, order, family, abbr
head -n 1 ABBR.csv > GENUS.csv
cat GENUS.tmp |
    sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> GENUS.csv

# clean
rm *.tmp *.bak

```

# Prepare sequences for lastz

```bash
cd ~/data/organelle/plastid/GENOMES

find . -maxdepth 1 -type d -path "*/*" |
    sort |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo >&2 "==> {}"
        
        if [ -e {}/chr.fasta ]; then
            echo >&2 "    {} has been processed"
            exit;
        fi

        egaz prepseq \
            {} \
            --gi -v --repeatmasker " --gff --parallel 8"
    '

# restore to original states
#for suffix in .2bit .fasta .fasta.fai .sizes .rm.out .rm.gff; do
#    find . -name "*${suffix}" | parallel --no-run-if-empty rm 
#done

```

# Aligning without outgroups

## Create alignments plans without outgroups

```bash
cd ~/data/organelle/plastid/summary

# tab-separated
# name  t   qs
cat GENUS.csv |
    grep -v "^#" |
    perl -na -F"," -e '
        BEGIN{
            $name = q{};
            %id_of = ();
        }

        chomp for @F;
        $F[4] =~ s/\W+/_/g;
        if ($F[4] ne $name) {
            if ($name) {
                my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
                my $t = shift @s;
                my $qs = join(q{,}, @s);
                printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
            }
            $name = $F[4];
            %id_of = ();
        }
        $id_of{$F[9]} = $F[0];

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = shift @s;
            my $qs = join(q{,}, @s);
            printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
        }
    ' \
    > genus.tsv

cat ABBR.csv |
    grep -v "^#" |
    perl -na -F"," -e '
        BEGIN{
            $name = q{};
            %id_of = ();
        }

        chomp for @F;
        $F[5] =~ s/\W+/_/g;
        if ($F[5] ne $name) {
            if ($name) {
                my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
                my $t = shift @s;
                my $qs = join(q{,}, @s);
                printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
            }
            $name = $F[5];
            %id_of = ();
        }
        $id_of{$F[9]} = $F[0]; # multiple chromosomes collapsed here

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = shift @s;
            my $qs = join(q{,}, @s);
            printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
        }
    ' \
    > family.tsv

```

```bash
cd ~/data/organelle/plastid/summary

cat <<'EOF' > egaz_template_multi.tt

# [% name %]
egaz template \
    ~/data/organelle/plastid/GENOMES/[% t %] \
[% FOREACH q IN qs -%]
    ~/data/organelle/plastid/GENOMES/[% q %] \
[% END -%]
[% IF o -%]
    ~/data/organelle/plastid/GENOMES/[% o %] \
    --outgroup [% o %] \
[% END -%]
    --multi -o [% name %] \
    --taxon ~/data/organelle/plastid/GENOMES/taxon_ncbi.csv \
    --rawphylo --aligndb --parallel 8 -v

EOF

# every genera
echo "mkdir -p ~/data/organelle/plastid/genus"  > ../cmd.txt
echo "cd       ~/data/organelle/plastid/genus" >> ../cmd.txt
cat genus.tsv |
    TT_FILE=egaz_template_multi.tt perl -MTemplate -nla -F"\t" -e '
        next unless scalar @F >= 3;
        
        my $tt = Template->new;
        $tt->process(
            $ENV{TT_FILE},
            {
                name       => $F[0],
                t          => $F[1],
                qs         => [ split /,/, $F[2] ],
                o          => $F[3],
            },
            \*STDOUT
        ) or die Template->error;

    ' \
    >> ../cmd.txt

# this is for finding outgroups
echo "mkdir -p ~/data/organelle/plastid/family"  > ../family.cmd.txt
echo "cd       ~/data/organelle/plastid/family" >> ../family.cmd.txt
cat family.tsv |
    TT_FILE=egaz_template_multi.tt perl -MTemplate -nla -F"\t" -e '
        next unless scalar @F >= 3;
        
        my $tt = Template->new;
        $tt->process(
            $ENV{TT_FILE},
            {
                name       => $F[0],
                t          => $F[1],
                qs         => [ split /,/, $F[2] ],
                o          => $F[3],
            },
            \*STDOUT
        ) or die Template->error;

    ' \
    >> ../family.cmd.txt

```

## Batch running for genera

```bash
mkdir -p ~/data/organelle/plastid/genus
cd ~/data/organelle/plastid/genus

bash ../cmd.txt 2>&1 | tee log_cmd.txt

#----------------------------#
# Step by step
#----------------------------#
# 1_pair
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 1_pair.sh | sort`; do
    echo "bash $f"
    echo
done > run_1.sh

# 2_rawphylo
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 2_rawphylo.sh | sort`; do
    echo "bash $f"
    echo
done > run_2.sh

# 3_multi
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 3_multi.sh | sort`; do
    echo "bash $f"
    echo
done > run_3.sh

# 6_chr_length
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 6_chr_length.sh | sort`; do
    echo "bash $f"
    echo
done > run_6.sh

# 7_multi_aligndb
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 7_multi_aligndb.sh | sort`; do
    echo "bash $f"
    echo
done > run_7.sh

cat run_1.sh | grep . | parallel -r -j 4  2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel -r -j 3  2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel -r -j 3  2>&1 | tee log_3.txt
cat run_6.sh | grep . | parallel -r -j 12 2>&1 | tee log_6.txt
cat run_7.sh | grep . | parallel -r -j 8  2>&1 | tee log_7.txt

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

## Alignments of families for outgroups.

```bash
mkdir -p ~/data/organelle/plastid/family
cd ~/data/organelle/plastid/family

time bash ../family.cmd.txt 2>&1 | tee log_cmd.txt

#----------------------------#
# Step by step
#----------------------------#
# 1_pair
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 1_pair.sh | sort`; do
    echo "bash $f"
    echo
done > run_1.sh

# 2_rawphylo
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 2_rawphylo.sh | sort`; do
    echo "bash $f"
    echo
done > run_2.sh

# 3_multi
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 3_multi.sh | sort`; do
    echo "bash $f"
    echo
done > run_3.sh

cat run_1.sh | grep . | parallel -r -j 4  2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel -r -j 3  2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel -r -j 3  2>&1 | tee log_3.txt

find ~/data/organelle/plastid/family -type f -name "*.nwk"

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw"   | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

# Aligning with outgroups

## Create `plastid_OG.md` for picking outgroups

Manually edit it then move to `~/Scripts/withncbi/doc/plastid_OG.md`.

```bash
cd ~/data/organelle/plastid/summary

cat GENUS.csv |
    grep -v "^#" |
    perl -na -F"," -e '
        BEGIN{
            ($phylum, $family, $genus, ) = (q{}, q{}, q{});
        }

        chomp for @F;

        if ($F[8] ne $phylum) {
            $phylum = $F[8];
            printf qq{\n# %s\n}, $phylum;
        }
        if ($F[5] ne $family) {
            $family = $F[5];
            printf qq{## %s\n}, $family;
        }
        $F[4] =~ s/\W+/_/g;
        if ($F[4] ne $genus) {
            $genus = $F[4];
            printf qq{%s\n}, $genus;
        }
    ' \
    > plastid_OG.md

```

## Create alignments plans with outgroups

```bash
cd ~/data/organelle/plastid/summary

# name  t   qs  o
cat genus.tsv |
    perl -nla -F"\t" -MPath::Tiny -e '
        BEGIN{
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{~/Scripts/withncbi/doc/plastid_OG.md})->lines({ chomp => 1});
            for (@ls) {
                @fs = split(/,/);
                $h{$fs[0]}= $fs[1];
            }
        }

        if (exists $h{$F[0]}) {
            printf qq{%s\t%s\t%s\t%s\n}, $F[0] . q{_OG}, $F[1], $F[2], $h{$F[0]};
        }
    ' \
    > genus_OG.tsv

# genera with outgroups
echo "mkdir -p ~/data/organelle/plastid/OG"  > ../OG.cmd.txt
echo "cd       ~/data/organelle/plastid/OG" >> ../OG.cmd.txt
cat genus_OG.tsv |
    TT_FILE=egaz_template_multi.tt perl -MTemplate -nla -F"\t" -e '
        next unless scalar @F >= 3;
        
        my $tt = Template->new;
        $tt->process(
            $ENV{TT_FILE},
            {
                name       => $F[0],
                t          => $F[1],
                qs         => [ split /,/, $F[2] ],
                o          => $F[3],
            },
            \*STDOUT
        ) or die Template->error;

    ' \
    >> ../OG.cmd.txt

```

In previous steps, we have manually edited `~/Scripts/withncbi/doc/plastid_OG.md` and generated
`genus_OG.tsv`.

*D* between target and outgroup should be around **0.05**.

```bash
mkdir -p ~/data/organelle/plastid/OG
cd ~/data/organelle/plastid/OG

time bash ../OG.cmd.txt 2>&1 | tee log_cmd.txt

#----------------------------#
# Step by step
#----------------------------#
# 1_pair
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 1_pair.sh | sort`; do
    echo "bash $f"
    echo
done > run_1.sh

# 2_rawphylo
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 2_rawphylo.sh | sort`; do
    echo "bash $f"
    echo
done > run_2.sh

# 3_multi
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 3_multi.sh | sort`; do
    echo "bash $f"
    echo
done > run_3.sh

# 6_chr_length
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 6_chr_length.sh | sort`; do
    echo "bash $f"
    echo
done > run_6.sh

# 7_multi_aligndb
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 7_multi_aligndb.sh | sort`; do
    echo "bash $f"
    echo
done > run_7.sh

cat run_1.sh | grep . | parallel -r -j 4  2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel -r -j 3  2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel -r -j 3  2>&1 | tee log_3.txt
cat run_6.sh | grep . | parallel -r -j 12 2>&1 | tee log_6.txt
cat run_7.sh | grep . | parallel -r -j 8  2>&1 | tee log_7.txt

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw"   | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

# Self alignments

```bash
cd ~/data/organelle/plastid/summary

cat <<'EOF' > egaz_templates_self.tt

# [% name %]
egaz template \
    ~/data/organelle/plastid/GENOMES/[% t %] \
[% FOREACH q IN qs -%]
    ~/data/organelle/plastid/GENOMES/[% q %] \
[% END -%]
    --self -o [% name %] \
    --taxon ~/data/organelle/plastid/GENOMES/taxon_ncbi.csv \
    --circos --parallel 8 -v

EOF

# every genera
echo "mkdir -p ~/data/organelle/plastid/self"  > ../self.cmd.txt
echo "cd       ~/data/organelle/plastid/self" >> ../self.cmd.txt
cat genus.tsv |
    TT_FILE=egaz_templates_self.tt perl -MTemplate -nla -F"\t" -e '
        next unless scalar @F >= 3;
        
        my $tt = Template->new;
        $tt->process(
            $ENV{TT_FILE},
            {
                name       => $F[0],
                t          => $F[1],
                qs         => [ split /,/, $F[2] ],
            },
            \*STDOUT
        ) or die Template->error;

    ' \
    >> ../self.cmd.txt

```

```bash
mkdir -p ~/data/organelle/plastid/self
cd ~/data/organelle/plastid/self

time bash ../self.cmd.txt 2>&1 | tee log_cmd.txt

#----------------------------#
# Step by step
#----------------------------#
# 1_self
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 1_self.sh | sort`; do
    echo "bash $f"
    echo
done > run_1.sh

# 3_proc
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 3_proc.sh | sort`; do
    echo "bash $f"
    echo
done > run_2.sh

# 4_circos
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 4_circos.sh | sort`; do
    echo "bash $f"
    echo
done > run_3.sh

cat run_1.sh | grep . | parallel -r -j 4  2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel -r -j 4  2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel -r -j 4  2>&1 | tee log_3.txt

# clean mysql
#find  /usr/local/var/mysql -type d -name "[A-Z]*" | parallel -r rm -fr

```

# LSC and SSC

IRA and IRB are presented by `plastid/self/${GENUS}/Results/${STRAIN}/${STRAIN}.links.tsv`.

```bash
find ~/data/organelle/plastid/self -type f -name "*.links.tsv" |
    xargs wc -l |
    sort -n |
    grep -v "total" |
    perl -nl -e 's/^\s*//g; /^(\d+)\s/ and print $1' |
    uniq -c |
    perl -nl -e '/^\s+(\d+)\s+(.+)$/ and print qq{$1\t$2}' |
    (echo -e 'count\tlines' && cat) |
    mlr --itsv --omd cat

```

Manually check strains not containing singular link.

| count | lines |
|:------|:------|
| 227   | 0     |
| 1949  | 1     |
| 44    | 2     |
| 22    | 3     |
| 4     | 4     |
| 2     | 10    |

| count | lines |
|------:|------:|
|   155 |     0 |
|  1709 |     1 |
|    36 |     2 |
|    17 |     3 |
|     1 |     4 |
|     1 |     5 |
|     1 |     8 |
|     1 |    10 |

There are 2 special strains (Asa_sieboldii, Epipo_aphyllum) which have no IR but a palindromic
sequence.

Create `ir_lsc_ssc.tsv` for slicing alignments.

Manually edit it then move to `~/Scripts/withncbi/doc/ir_lsc_ssc.tsv`.

`#genus abbr role accession chr_size IR LSC SSC`

* `NA` - not self-aligned
* `MULTI` - multiply link records
* `NONE` - no link records
* `WRONG` - unexpected link records

```bash
cd ~/data/organelle/plastid/summary

cat ABBR.csv |
    grep -v "^#" |
    perl -nla -F"," -MAlignDB::IntSpan -MPath::Tiny -e '
        BEGIN {
            %seen = ();
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{family.tsv})->lines({ chomp => 1});
            for (@ls) {
                $seen{$_}++ for split(/,|\s+/);
            }

            %target = ();
            %queries = ();
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{genus.tsv})->lines({ chomp => 1});
            for (@ls) {
                $target{(split /\t/)[1]}++;
                $queries{$_}++ for split(/,/, (split /\t/)[2]);
            }

            %outgroup = ();
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{genus_OG.tsv})->lines({ chomp => 1});
            for (@ls) {
                $outgroup{(split /\t/)[3]}++;
            }
        }

        chomp for @F;

        my $genus = $F[4];
        my $abbr = $F[9];

        next unless $seen{$abbr};

        my $role = "UNKNOWN";
        if ($target{$abbr}) {
            $role = "Target";
        }
        elsif ($queries{$abbr}) {
            $role = "Queries";
        }
        elsif ($outgroup{$abbr}) {
            $role = "Outgroup";
        }

        my $size_file = qq{$ENV{HOME}/data/organelle/plastid_self.working/$genus/Genomes/$abbr/chr.sizes};
        my $link_file = qq{$ENV{HOME}/data/organelle/plastid_self.working/$genus/Results/$abbr/$abbr.links.tsv};

        if (! -e $size_file or ! -e $link_file) {
            if ($outgroup{$abbr}) {
                print q{#} . join(qq{\t}, $genus, $abbr, $role, (q{NA}) x 5 );
            }
            next;
        }

        my $accession = `cat $size_file | cut -f 1 | xargs echo`;
        my $chr_size = `cat $size_file | cut -f 2 | xargs echo`;
        my $link = `cat $link_file | xargs echo`;
        chomp $accession; chomp $chr_size; chomp $link;

        if (split(q{ }, $chr_size) > 2 or split(q{ }, $link) > 2) {
            print q{#} . join(qq{\t}, $genus, $abbr, $role, $accession, $chr_size, (q{MULTI}) x 3 );
            next;
        }

        if ($link !~ m{:(\d+\-\d+)\s+.+:(\d+\-\d+)$}) {
            print q{#} . join(qq{\t}, $genus, $abbr, $role, $accession, $chr_size, (q{NONE}) x 3 );
            next;
        }

        my $ira = AlignDB::IntSpan->new($1);
        my $irb = AlignDB::IntSpan->new($2);
        my $ir = AlignDB::IntSpan->new->add($ira)->add($irb);

        if ($ira->trim(1)->contains(1)
            or $irb->trim(1)->contains(1)
            or $ira->trim(1)->contains($chr_size)
            or $irb->trim(1)->contains($chr_size)
            or ($ira->max + 1 > $irb->min - 1)) {
            print q{#} . join(qq{\t}, $genus, $abbr, $role, $accession, $chr_size, (q{WRONG}) x 3 );
            next;
        }

        my $chr = AlignDB::IntSpan->new->add_pair(1, $chr_size);
        my $interval = AlignDB::IntSpan->new->add_pair($ira->max + 1, $irb->min - 1);

        my $d_s_a = $ira->distance(AlignDB::IntSpan->new(1));
        my $d_b_e = $irb->distance(AlignDB::IntSpan->new($chr_size));
        my $d_a_b = $ira->distance($irb);
        $d_s_a = 0 if $d_s_a < 0;
        $d_b_e = 0 if $d_b_e < 0;

        my $lsc = AlignDB::IntSpan->new;
        my $ssc = AlignDB::IntSpan->new;

        if ($d_s_a + $d_b_e > $d_a_b) { # LSC
            $ssc = $interval->copy;
            $lsc = $chr->diff($ir)->diff($ssc);
        }
        else  { # SSC
            $lsc = $interval->copy;
            $ssc = $chr->diff($ir)->diff($lsc);
        }

        print join(qq{\t}, $genus, $abbr, $role, $accession, $chr_size, $ir, $lsc, $ssc );
    ' \
    > ir_lsc_ssc.tsv

```

## Can't get clear IR information

* Grateloupia
  * Grat_filicina
  * Grat_taiwanensis
* Caulerpa
  * Cau_cliftonii
  * Cau_racemosa
* Caloglossa
  * Calog_beccarii
  * Calog_intermedia
  * Calog_monosticha
* Pisum
  * Pisum_fulvum
  * Pisum_sativum
  * Pisum_sativum_subsp_elatius
* Dasya
  * Dasya_naccarioides
* Diplazium
  * Diplazium_unilobum
* Bryopsis
  * Bryop_plumosa
  * Bryop_sp_HV04063
  * Bry_sp_HV04063
* Medicago
  * Med_falcata
  * Med_hybrida
  * Med_papillosa
  * Med_truncatula
* Aegilops
  * Aeg_cylindrica
  * Aeg_geniculata
  * Aeg_speltoides
  * Aeg_tauschii
* Prototheca
  * Prot_cutis
  * Prot_stagnorum
  * Prot_zopfii
* Cryptomonas
  * Cryptomo_curvata
  * Cryptomo_paramecium
* Monotropa
  * Monotropa_hypopitys
* Liagora
  * Liagora_brachyclada
  * Liagora_harveyana
* Taiwania
  * Tai_cryptomerioides
  * Tai_flousiana
* Pinus
  * Pinus_armandii
  * Pinus_bungeana
  * Pinus_contorta
  * Pinus_gerardiana
  * Pinus_greggii
  * Pinus_jaliscana
  * Pinus_koraiensis
  * Pinus_krempfii
  * Pinus_lambertiana
  * Pinus_longaeva
  * Pinus_massoniana
  * Pinus_monophylla
  * Pinus_nelsonii
  * Pinus_oocarpa
  * Pinus_pinea
  * Pinus_sibirica
  * Pinus_strobus
  * Pinus_sylvestris
  * Pinus_tabuliformis
  * Pinus_taeda
  * Pinus_taiwanensis
  * Pinus_thunbergii
* Taxus
  * Taxus_mairei
* Picea
  * Pic_abies
  * Pic_asperata
  * Pic_crassifolia
  * Pic_glauca
  * Pic_jezoensis
  * Pic_morrisonicola
  * Pic_sitchensis
* Gracilariopsis
  * Gracilario_chorda
  * Gracilario_lemaneiformis
* Fragaria
  * Frag_mandshurica
  * Frag_vesca_subsp_bracteata
* Ulva
  * Ulva_fasciata
  * Ulva_flexuosa
  * Ulva_linza
  * Ulva_prolifera
* Monomorphina
  * Monom_parapyrum
* Epipogium
  * Epipo_aphyllum
* Euglena
  * Euglena_archaeoplastidiata
  * Euglena_viridis
* Chlorella
  * Chlore_heliozoae
  * Chlore_sorokiniana
  * Chlore_variabilis
  * Chlore_vulgaris
* Erodium
  * Ero_carvifolium
  * Ero_crassifolium
  * Ero_manescavi
  * Ero_rupestre
* Larix
  * Lar_decidua
  * Lar_sibirica
* Amentotaxus
  * Ame_argotaenia
  * Ame_formosana
* Pyropia
  * Pyro_perforata
* Ceramium
  * Ceram_cimbricum
  * Ceram_japonicum
  * Ceram_sungminbooi
* Hildenbrandia
  * Hilde_rivularis
  * Hilde_rubra
* Pilostyles
  * Pilo_aethiopica
  * Pilo_hamiltonii
* Codium
  * Codi_decorticatum
  * Codi_sp_arenicola
* Torreya
  * Torreya_fargesii
  * Torreya_grandis
* Vertebrata
  * Vert_australis
  * Vert_isogona
  * Vert_lanosa
  * Vert_thuyoides
* Wisteria
  * Wis_floribunda
  * Wis_sinensis
* Phelipanche
  * Pheli_purpurea
  * Pheli_ramosa
* Glycyrrhiza
  * Glycy_glabra
  * Glycy_lepidota
* Cephalotaxus
  * Cephalo_wilsoniana
* Polysiphonia
  * Polysi_brodiei
  * Polysi_elongata
  * Polysi_infestans
  * Polysi_schneideri
  * Polysi_scopulorum
  * Polysi_sertularioides
  * Polysi_stricta
* Lathyrus
  * Lathy_clymenum
  * Lathy_davidii
  * Lathy_graminifolius
  * Lathy_inconspicuus
  * Lathy_littoralis
  * Lathy_ochroleucus
  * Lathy_odoratus
  * Lathy_palustris
  * Lathy_sativus
  * Lathy_tingitanus
* Babesia
  * Bab_orientalis
* Gelidium
  * Gelidi_elegans
  * Gelidi_vagum
* Bostrychia
  * Bos_moritziana
  * Bos_simpliciuscula
  * Bos_tenella
* Membranoptera
  * Mem_platyphylla
  * Mem_tenuis
  * Mem_weeksiae
* Astragalus
  * Astra_mongholicus
  * Astra_mongholicus_var_nakaianus
* Asarum
  * Asa_minus
  * Asa_sieboldii
* Gracilaria
  * Gracilaria_changii
  * Gracilaria_chilensis
  * Gracilaria_firma
  * Gracilaria_salicornia
  * Gracilaria_tenuistipitata_var_liui
  * Gracilaria_vermiculophylla
* Plasmodium
  * Plas_chabaudi_chabaudi
  * Plas_falciparum_HB3
  * Plas_gallinaceum
  * Plas_relictum
  * Plas_vivax
* Triticum
  * Trit_urartu
* Trifolium
  * Trif_aureum
  * Trif_boissieri
  * Trif_glanduliferum
  * Trif_grandiflorum
  * Trif_strictum

There are 2 special strains which has only one palindromic sequence rather than IR. (as mentioned
before)

## Slices of IR, LSC and SSC

Without outgroups.

Be cautious to alignments with low coverage.

```bash
mkdir -p ~/data/organelle/plastid/slices
cd ~/data/organelle/plastid/slices

cat ~/Scripts/withncbi/doc/ir_lsc_ssc.tsv |
    perl -nla -F"\t" -MAlignDB::IntSpan -Mstrict -Mwarnings -e '
        /^#/ and next;
        $F[2] eq q{Target} or next;
        $F[5] =~ /\d+/ or next;

        print qq{# $F[0]};

        next unless -e "$ENV{HOME}/data/organelle/plastid/genus/$F[0]/$F[0]_refined/$F[3].synNet.maf.gz.fas.gz";

        my %rl_of = (
            IR  => $F[5],
            LSC => $F[6],
            SSC => $F[7],
        );

        my $lsc = AlignDB::IntSpan->new($F[6]);
        # next unless $lsc->span_size == 2; # for testing

        my $max_seg = 3;
        my $segment = int($lsc->size / $max_seg /2);
        for my $i (1 .. $max_seg) {
            my $slice_5 = AlignDB::IntSpan->new;
            my $slice_3 = AlignDB::IntSpan->new;

            # these are all indexes of LCS
            my $seg_5_start = $segment * ($i - 1) + 1;
            my $seg_5_end = $segment * $i;
            my $seg_3_start = $lsc->size - $segment * $i + 1;
            my $seg_3_end = $lsc->size - $segment * ($i - 1);

            if ($lsc->span_size == 1) {
                $slice_5 = $lsc->slice($seg_5_start, $seg_5_end);
                $slice_3 = $lsc->slice($seg_3_start, $seg_3_end);
            }
            elsif ($lsc->span_size == 2) { # start point inside LSC
                my ($lsc1, $lsc2) = $lsc->sets;

                if ($lsc2->size >= $seg_5_end) {
                    $slice_5 = $lsc2->slice($seg_5_start, $seg_5_end);
                }
                elsif ($lsc2->size >= $seg_5_start) {
                    $slice_5 = $lsc2->slice($seg_5_start, $lsc2->size);
                    $slice_5->add( $lsc1->slice(1, $segment - ($lsc2->size - $seg_5_start + 1) ) );
                }
                else {
                    $slice_5 = $lsc1->slice($seg_5_start - $lsc2->size, $seg_5_end - $lsc2->size);
                }

                if ($lsc1->size >= $lsc->size - $seg_3_start) {
                    $slice_3 = $lsc1->slice($seg_3_start - $lsc2->size, $seg_3_end - $lsc2->size);
                }
                elsif ($lsc1->size >= $lsc->size - $seg_3_end) {
                    $slice_3 = $lsc1->slice(1, $seg_3_end - $lsc2->size);
                    $slice_3->add( $lsc2->slice( $seg_3_start - ($seg_3_end - $lsc2->size), $lsc2->size ) );
                }
                else {
                    $slice_3 = $lsc2->slice($seg_3_start, $seg_3_end);
                }
            }

            my $slice = AlignDB::IntSpan->new;
            $slice->add($slice_5)->add($slice_3);
            $slice = $slice->fill(10); # fill small gaps in LSC3
            $rl_of{"LSC$i"} = $slice->as_string;
        }

        for my $key (sort keys %rl_of) {
            print qq{jrunlist cover <(echo $F[3]:$rl_of{$key}) -o $F[0].$key.yml};
            print qq{fasops slice -n $F[1] -o $F[0].$key.fas \\};
            print qq{    ~/data/organelle/plastid/genus/$F[0]/$F[0]_refined/$F[3].synNet.maf.gz.fas.gz \\};
            print qq{    $F[0].$key.yml};
            print qq{perl ~/Scripts/alignDB/alignDB.pl \\};
            print qq{    -d $F[0]_$key \\};
            print qq{    -da ~/data/organelle/plastid_slices/$F[0].$key.fas \\};
            print qq{    -a ~/data/organelle/plastid/genus/$F[0]/Stats/anno.yml \\};
            print qq{    -chr ~/data/organelle/plastid/genus/$F[0]/chr_length.csv \\};
            print qq{    --lt 1000 --parallel 8 --batch 5 \\};
            print qq{    --run common};
            print qq{};
        }

        print qq{};
    ' \
    > slices.sh

```

Run the generated bash file.

```bash
cd ~/data/organelle/plastid_slices

bash slices.sh
perl ~/Scripts/fig_table/collect_common_basic.pl -d .

```

# Cyanobacteria

## Genus and Species counts

```sql
# Genus
SELECT
    genus_id, genus, COUNT(*) strain_count
FROM
    gr_prok.gr
WHERE
    `group` = 'Cyanobacteria'
        AND status NOT IN ('Contig' , 'Scaffold')
GROUP BY genus_id
HAVING strain_count > 1
ORDER BY genus
```

| genus_id | genus               | strain_count |
|---------:|:--------------------|-------------:|
|     1163 | Anabaena            |            3 |
|    35823 | Arthrospira         |            2 |
|     1186 | Calothrix           |            3 |
|   102234 | Cyanobacterium      |            2 |
|    43988 | Cyanothece          |            6 |
|    33071 | Gloeobacter         |            2 |
|     1125 | Microcystis         |            2 |
|     1177 | Nostoc              |            4 |
|     1158 | Oscillatoria        |            2 |
|     1218 | Prochlorococcus     |           14 |
|     1129 | Synechococcus       |           20 |
|     1142 | Synechocystis       |            5 |
|   146785 | Thermosynechococcus |            2 |

```sql
# Species
SELECT
    species_id, species, COUNT(*) strain_count
FROM
    gr_prok.gr
WHERE
    `group` = 'Cyanobacteria'
        AND status NOT IN ('Contig' , 'Scaffold')
GROUP BY species_id
HAVING strain_count > 1
ORDER BY species
```

| species_id | species                    | strain_count |
|-----------:|:---------------------------|-------------:|
|     118562 | Arthrospira platensis      |            2 |
|       1126 | Microcystis aeruginosa     |            2 |
|       1219 | Prochlorococcus marinus    |           12 |
|      32046 | Synechococcus elongatus    |            2 |
|       1148 | Synechocystis sp. PCC 6803 |            4 |

# Summary

## Copy xlsx files

```bash
mkdir -p ~/data/organelle/plastid/summary/xlsx
cd ~/data/organelle/plastid/summary/xlsx

find ../../genus -type f -name "*.common.xlsx" |
    grep -v "vs[A-Z]" |
    parallel cp {} .

find ../../OG -type f -name "*.common.xlsx" |
    grep -v "vs[A-Z]" |
    parallel cp {} .

```

## Genome list

Create `list.csv` from `GENUS.csv` with sequence lengths.

```bash
mkdir -p ~/data/organelle/plastid/summary/table
cd ~/data/organelle/plastid/summary/table

# manually set orders in `plastid_OG.md`
perl -l -MPath::Tiny -e '
    BEGIN {
        @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_} 
            map {s/,\w+//; $_} 
            map {s/^###\s*//; $_} 
            path(q{~/Scripts/withncbi/doc/plastid_OG.md})->lines( { chomp => 1}); 
    }
    for (@ls) { 
        (/^\s*$/ or /^##\s+/ or /^#\s+(\w+)/) and next; 
        print $_
    }
    ' \
    > genus_all.lst

# abbr accession length
find ../../genus -type f -name "chr_length.csv" |
    parallel --jobs 1 --keep-order -r '
        perl -nl -e '\''
            BEGIN {
                our $l = { }; # avoid parallel replace string
            }
            
            next unless /\w+,\d+/;
            my ($common_name, undef, $chr, $length) = split /,/;
            if (exists $l->{$common_name}) {
                $l->{$common_name}{$chr} = $length;
            }
            else {
                $l->{$common_name} = {$chr => $length};
            }
            
            END {
                for my $common_name (keys %{$l}) {
                    my $chrs = join "|", sort keys %{$l->{$common_name}};
                    my $length = 0;
                    $length += $_ for values %{$l->{$common_name}};
                    print qq{$common_name\t$chrs\t$length}
                }
            }
        '\'' {}
    ' \
    > length.tmp
cat length.tmp | datamash check
#2248 lines, 3 fields

# phylum family genus abbr taxon_id
cat ~/data/organelle/plastid/summary/GENUS.csv |
    grep -v "^#" |
    perl -nla -F"," -e 'print join qq{\t}, ($F[8], $F[5], $F[4], $F[9], $F[0], )' |
    sort |
    uniq \
    > abbr.tmp
cat abbr.tmp | datamash check
#2250 lines, 5 fields

tsv-join \
    abbr.tmp \
    --data-fields 4 \
    -f length.tmp \
    --key-fields 1 \
    --append-fields 2,3 \
    > list.tmp
cat list.tmp | datamash check
#2248 lines, 7 fields

# sort as orders in plastid_OG.md
echo -e "#phylum,family,genus,abbr,taxon_id,accession,length" > list.csv
cat list.tmp |
    perl -nl -a -MPath::Tiny -e '
        BEGIN{
            %genus;
            my @l1 = path(q{genus_all.lst})->lines({ chomp => 1});
            $genus{$l1[$_]} = $_ for (0 .. $#l1);
        }
        my $idx = $genus{$F[2]};
        die qq{$_\n} unless defined $idx;
        print qq{$_\t$idx};
    ' |
    sort -n -k8,8 |
    cut -f 1-7 |
    tr $'\t' ',' \
    >> list.csv
cat list.csv | datamash check -t,
#2249 lines, 7 fields

rm *.tmp

```

## Statistics of genome alignments

Some genera will be filtered out here.

Criteria:

* Coverage >= 0.5
* Total number of indels >= 100
* Genome D < 0.05

```bash
cd ~/data/organelle/plastid/summary/xlsx

cat <<'EOF' > Table_alignment.tt
---
autofit: A:F
texts:
  - text: "Genus"
    pos: A1
  - text: "No. of genomes"
    pos: B1
  - text: "Aligned length (Mb)"
    pos: C1
  - text: "Indels Per 100 bp"
    pos: D1
  - text: "D on average"
    pos: E1
  - text: "GC-content"
    pos: F1
[% FOREACH item IN data -%]
  - text: [% item.name %]
    pos: A[% loop.index + 2 %]
[% END -%]
borders:
  - range: A1:F1
    top: 1
  - range: A1:F1
    bottom: 1
ranges:
[% FOREACH item IN data -%]
  [% item.file %]:
    basic:
      - copy: B2
        paste: B[% loop.index + 2 %]
      - copy: B4
        paste: C[% loop.index + 2 %]
      - copy: B5
        paste: D[% loop.index + 2 %]
      - copy: B7
        paste: E[% loop.index + 2 %]
      - copy: B8
        paste: F[% loop.index + 2 %]
[% END -%]
EOF

cat ~/data/organelle/plastid/summary/table/genus_all.lst |
    grep -v "^#" |
    TT_FILE=Table_alignment.tt perl -MTemplate -nl -e '
        push @data, { name => $_, file => qq{$_.common.xlsx}, }; 
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, })
                or die Template->error;
        }
    ' \
    > Table_alignment_all.yml

perl ~/Scripts/fig_table/xlsx_table.pl -i Table_alignment_all.yml
perl ~/Scripts/fig_table/xlsx2csv.pl -f Table_alignment_all.xlsx > Table_alignment_all.csv

cp -f Table_alignment_all.xlsx ~/data/organelle/plastid/summary/table
cp -f Table_alignment_all.csv ~/data/organelle/plastid/summary/table

```

```bash
cd ~/data/organelle/plastid/summary/table

echo "Genus,avg_size" > group_avg_size.csv
cat list.csv |
    grep -v "#" |
    perl -nla -F, -e '
        $count{$F[2]}++;
        $sum{$F[2]} += $F[6];
        END {
            for $k (sort keys %count) {
                printf qq{%s,%d\n}, $k, $sum{$k}/$count{$k};
            }
        }
    ' \
    >> group_avg_size.csv

cat Table_alignment_all.csv group_avg_size.csv |
    perl ~/Scripts/withncbi/util/merge_csv.pl -f 0 --concat -o stdout \
    > Table_alignment_all.1.csv

echo "Genus,coverage" > group_coverage.csv
cat Table_alignment_all.1.csv |
    perl -nla -F',' -e '
        $F[2] =~ /[\.\d]+/ or next;
        $F[6] =~ /[\.\d]+/ or next;
        $c = $F[2] * 1000 * 1000 / $F[6];
        print qq{$F[0],$c};
    ' \
    >> group_coverage.csv

cat Table_alignment_all.1.csv group_coverage.csv |
    perl ~/Scripts/withncbi/util/merge_csv.pl -f 0 --concat -o stdout \
    > Table_alignment_all.2.csv

echo "Genus,indels" > group_indels.csv
cat Table_alignment_all.2.csv |
    perl -nla -F',' -e '
        $F[6] =~ /[\.\d]+/ or next;
        $c = $F[3] / 100 * $F[2] * 1000 * 1000;
        print qq{$F[0],$c};
    ' \
    >> group_indels.csv

cat Table_alignment_all.2.csv group_indels.csv |
    perl ~/Scripts/withncbi/util/merge_csv.pl -f 0 --concat -o stdout |
    grep -v ',,' \
    > Table_alignment_for_filter.csv

# real filter
cat Table_alignment_for_filter.csv |
    perl -nla -F',' -e '
        $F[6] =~ /[\.\d]+/ or next;
        $F[0] =~ s/"//g;
        print $F[0] if ($F[7] >= 0.5 and $F[8] >= 100 and $F[4] < 0.05);
    ' \
    > genus.lst

rm ~/data/organelle/plastid/summary/table/Table_alignment_all.[0-9].csv
rm ~/data/organelle/plastid/summary/table/group_*csv

#
cd ~/data/organelle/plastid/summary/xlsx
cat ~/data/organelle/plastid/summary/table/genus.lst |
    grep -v "^#" |
    TT_FILE=Table_alignment.tt perl -MTemplate -nl -e '
        push @data, { name => $_, file => qq{$_.common.xlsx}, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, })
                or die Template->error;
        }
    ' \
    > Table_alignment.yml

perl ~/Scripts/fig_table/xlsx_table.pl -i Table_alignment.yml
perl ~/Scripts/fig_table/xlsx2csv.pl -f Table_alignment.xlsx > Table_alignment.csv

cp -f ~/data/organelle/plastid/summary/xlsx/Table_alignment.xlsx ~/data/organelle/plastid/summary/table
cp -f ~/data/organelle/plastid/summary/xlsx/Table_alignment.csv ~/data/organelle/plastid/summary/table

```

## Groups

```bash
mkdir -p ~/data/organelle/plastid/summary/group
cd ~/data/organelle/plastid/summary/group

perl -l -MPath::Tiny -e '
    BEGIN {
        @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_}
            map {s/,\w+//; $_}
            map {s/^###\s*//; $_}
            path(q{~/Scripts/withncbi/doc/plastid_OG.md})->lines( { chomp => 1});
    }
    for (@ls) {
        (/^\s*$/ or /^##\s+/) and next;
        if (/^#\s+(\w+)/) {
            $fh = path("$1.txt")->openw;
            next;
        } else {
            print {$fh} $_;
        }
    }
    '

grep -Fx -f ../table/genus.lst Angiosperms.txt > Angiosperms.lst

find . -type f -name "*.txt" |
    xargs cat |
    grep -v -Fx -f Angiosperms.txt |
    grep -Fx -f ../table/genus.lst \
    > Others.lst

cat ../table/Table_alignment.csv |
    cut -d, -f 1,5 |
    perl -nla -F',' -e '$F[1] > 0.02 and $F[1] <= 0.05 and print $F[0];' |
    grep -Fx -f Angiosperms.lst \
    > group_3.lst

cat ../table/Table_alignment.csv |
    cut -d, -f 1,5 |
    perl -nla -F',' -e '$F[1] > 0.005 and $F[1] <= 0.02 and print $F[0];' |
    grep -Fx -f Angiosperms.lst \
    > group_2.lst

cat ../table/Table_alignment.csv |
    cut -d, -f 1,5 |
    perl -nla -F',' -e '$F[1] <= 0.005 and print $F[0];' |
    grep -Fx -f Angiosperms.lst \
    > group_1.lst

rm *.txt

```

NCBI Taxonomy tree

```bash
mkdir -p ~/data/organelle/plastid/summary/group
cd ~/data/organelle/plastid/summary/group

cat Others.lst |
    grep -v "^#" |
    perl -e '
        @ls = <>;
        $str = qq{bp_taxonomy2tree.pl \\\n};
        for (@ls) {
            chomp;
            $str .= qq{    -s "$_" \\\n};
        }
        $str .= qq{    -e \n};
        print $str;
    ' \
    > Others.tree.sh

bash Others.tree.sh >Others.nwk

nw_display -w 600 -s Others.nwk |
    rsvg-convert -f pdf -o Others.pdf

```

## Phylogenic trees of each genus with outgroup

```bash
mkdir -p ~/data/organelle/plastid/summary/trees
cd ~/data/organelle/plastid/summary/trees

cat ~/Scripts/withncbi/doc/plastid_OG.md |
    grep -v "^#" |
    grep . |
    cut -d',' -f 1 \
    > list.txt

find ../../OG -type f -path "*Results*" -name "*.nwk" |
    grep -v ".raw." |
    parallel -j 1 cp {} trees

```

## d1, d2

`collect_xlsx.pl`

```bash
cd ~/data/organelle/plastid/summary/xlsx

cat <<'EOF' > cmd_collect_d1_d2.tt
perl ~/Scripts/fig_table/collect_xlsx.pl \
[% FOREACH item IN data -%]
    -f [% item.name %].common.xlsx \
    -s d1_pi_gc_cv \
    -n [% item.name %] \
[% END -%]
    -o cmd_d1.xlsx

perl ~/Scripts/fig_table/collect_xlsx.pl \
[% FOREACH item IN data -%]
    -f [% item.name %].common.xlsx \
    -s d2_pi_gc_cv \
    -n [% item.name %] \
[% END -%]
    -o cmd_d2.xlsx

perl ~/Scripts/fig_table/collect_xlsx.pl \
[% FOREACH item IN data -%]
    -f [% item.name %].common.xlsx \
    -s d1_comb_pi_gc_cv \
    -n [% item.name %] \
[% END -%] 
    -o cmd_d1_comb.xlsx

perl ~/Scripts/fig_table/collect_xlsx.pl \
[% FOREACH item IN data -%]
    -f [% item.name %].common.xlsx \
    -s d2_comb_pi_gc_cv \
    -n [% item.name %] \
[% END -%] 
    -o cmd_d2_comb.xlsx

EOF

cat ../table/genus.lst |
    grep -v "^#" |
    TT_FILE=cmd_collect_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, }) 
                or die Template->error;
        }
    ' \
    > cmd_collect_d1_d2.sh

bash cmd_collect_d1_d2.sh

```

`sep_chart.pl`

```bash
mkdir -p ~/data/organelle/plastid/summary/fig

cd ~/data/organelle/plastid/summary/xlsx

cat <<'EOF' > cmd_chart_d1_d2.tt
perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d1.xlsx \
    -xl "Distance to indels ({italic(d)[1]})" \
    -yl "Nucleotide divergence ({italic(D)})" \
    -xr "A2:A8" -yr "B2:B8" \
    --y_min 0.0 --y_max [% y_max %] \
    -x_min 0 -x_max 5 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot -ms

perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d1_comb.xlsx \
    -xl "Distance to indels ({italic(d)[1]})" \
    -yl "Nucleotide divergence ({italic(D)})" \
    -xr "A2:A8" -yr "B2:B8" \
    --y_min 0.0 --y_max [% y_max %] \
    -x_min 0 -x_max 5 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot

perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d2.xlsx \
    -xl "Reciprocal of indel density ({italic(d)[2]})" \
    -yl "Nucleotide divergence ({italic(D)})" \
    -xr "A2:A23" -yr "B2:B23" \
    --y_min 0.0 --y_max [% y_max2 %] \
    -x_min 0 -x_max 20 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot -ms

perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d2_comb.xlsx \
    -xl "Reciprocal of indel density ({italic(d)[2]})" \
    -yl "Nucleotide divergence ({italic(D)})" \
    -xr "A2:A23" -yr "B2:B23" \
    --y_min 0.0 --y_max [% y_max2 %] \
    -x_min 0 -x_max 20 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot

EOF

cat ../group/group_1.lst |
    TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.01,
                y_max2 => 0.01,
                postfix => q{group_1}, }) 
                or die Template->error;
        }
    ' \
    > cmd_chart_group_1.sh

cat ../group/group_2.lst |
    TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.03,
                y_max2 => 0.03,
                postfix => q{group_2}, }) 
                or die Template->error;
        }
    ' \
    > cmd_chart_group_2.sh

cat ../group/group_3.lst |
    TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.05,
                y_max2 => 0.05,
                postfix => q{group_3}, }) 
                or die Template->error;
        }
    ' \
    > cmd_chart_group_3.sh

cat ../group/Others.lst |
    TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.15,
                y_max2 => 0.15,
                postfix => q{Others}, }) 
                or die Template->error;
        }
    ' \
    > cmd_chart_Others.sh

bash cmd_chart_group_1.sh
bash cmd_chart_group_2.sh
bash cmd_chart_group_3.sh
bash cmd_chart_Others.sh

rm ../xlsx/*.csv
mv ../xlsx/*.pdf ../fig

# Coreldraw doesn't play well with computer modern fonts (latex math).
# perl ~/Scripts/fig_table/tikz_chart.pl -i cmd_plastid_d1_A2A8_B2B8.group_1.csv -xl 'Distance to indels ($d_1$)' -yl 'Nucleotide divergence ($D$)' --y_min 0.0 --y_max 0.01 -x_min 0 -x_max 5 --style_dot --pdf

```

