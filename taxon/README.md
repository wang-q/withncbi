# Process plastid genomes

The following command lines are about how I processed the plastid
genomes of green plants. Many tools of `taxon/` are used here, which
makes a good example for users.

## Work flow.

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
Save page to a local file, html only. In this case, it's `doc/green_plants_plastid_181127.html`.

All [Eukaryota](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?opt=plastid&taxid=2759),
`doc/eukaryota_plastid_181127.html`.

Or [this link](http://www.ncbi.nlm.nih.gov/genome/browse/?report=5).

```text
Eukaryota (2759)                2938
    Viridiplantae (33090)       2727
        Chlorophyta (3041)      121
        Streptophyta (35493)    2606
```

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions from all history
pages.

```csv
id,acc
996148,NC_017006
```

Got **3358** accessions.

```bash
mkdir -p ~/data/organelle/plastid/GENOMES
cd ~/data/organelle/plastid/GENOMES

rm webpage_id_seq.csv
perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/eukaryota_plastid_181127.html \
    >> webpage_id_seq.csv    

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/green_plants_plastid_181127.html \
    >> webpage_id_seq.csv    

```

Use `taxon/gb_taxon_locus.pl` to extract information from refseq genbank files.

```bash
cd ~/data/organelle/plastid/GENOMES

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.1.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.2.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/plastid.3.genomic.gbff.gz

gzip -dcf plastid.*.genomic.gbff.gz > genomic.gbff

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv

rm genomic.gbff

# 3336
cat refseq_id_seq.csv | grep -v "^#" | wc -l

# combine
cat webpage_id_seq.csv refseq_id_seq.csv |
    sort -u | # duplicated id-seq pair
    sort -t, -k1,1 \
    > id_seq.csv

# 3358
cat id_seq.csv | grep -v "^#" | wc -l

```

# Add lineage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `CHECKME.csv`.

```bash
mkdir -p ~/data/organelle/plastid/summary
cd ~/data/organelle/plastid/summary

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
cd ~/data/organelle/plastid/summary

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
sed -i".bak" "s/Phaeophyceae,NA/Phaeophyceae,Ochrophyta/" CHECKME.csv
sed -i".bak" "s/Eustigmatophyceae,NA/Eustigmatophyceae,Ochrophyta/" CHECKME.csv
sed -i".bak" "s/Xanthophyceae,NA/Xanthophyceae,Ochrophyta/" CHECKME.csv
sed -i".bak" "s/Pelagophyceae,NA/Pelagophyceae,Ochrophyta/" CHECKME.csv
sed -i".bak" "s/Raphidophyceae,NA/Raphidophyceae,Ochrophyta/" CHECKME.csv
sed -i".bak" "s/Synurophyceae,NA/Synurophyceae,Ochrophyta/" CHECKME.csv
sed -i".bak" "s/NA,Phaeophyceae/Phaeophyceae,Ochrophyta/" CHECKME.csv

# Rhodophyta 红藻门
sed -i".bak" "s/Bangiophyceae,NA/Bangiophyceae,Rhodophyta/" CHECKME.csv
sed -i".bak" "s/Compsopogonophyceae,NA/Compsopogonophyceae,Rhodophyta/" CHECKME.csv
sed -i".bak" "s/Florideophyceae,NA/Florideophyceae,Rhodophyta/" CHECKME.csv
sed -i".bak" "s/Rhodellophyceae,NA/Rhodellophyceae,Rhodophyta/" CHECKME.csv
sed -i".bak" "s/Stylonematophyceae,NA/Stylonematophyceae,Rhodophyta/" CHECKME.csv

# Cryptophyta 隐藻门
sed -i".bak" "s/Cryptomonadales,NA,NA/Cryptomonadales,Cryptophyceae,Cryptophyta/" CHECKME.csv
sed -i".bak" "s/Pyrenomonadales,NA,NA/Pyrenomonadales,Cryptophyceae,Cryptophyta/" CHECKME.csv
sed -i".bak" "s/Cryptophyta,NA/Cryptophyta,Cryptophyta/" CHECKME.csv

# Charophyta 轮藻门
sed -i".bak" "s/Charophyceae,Streptophyta/Charophyceae,Charophyta/" CHECKME.csv
sed -i".bak" "s/Chlorokybophyceae,Streptophyta/Chlorokybophyceae,Charophyta/" CHECKME.csv
sed -i".bak" "s/Coleochaetophyceae,Streptophyta/Coleochaetophyceae,Charophyta/" CHECKME.csv
sed -i".bak" "s/Zygnemophyceae,Streptophyta/Zygnemophyceae,Charophyta/" CHECKME.csv

# Chlorophyta 绿藻门
sed -i".bak" "s/Mesostigmatophyceae,Streptophyta/Mesostigmatophyceae,Chlorophyta/" CHECKME.csv

# Bryophytes
sed -i".bak" "s/Marchantiopsida,Streptophyta/Marchantiopsida,Bryophytes/" CHECKME.csv
sed -i".bak" "s/Leiosporocerotopsida,Streptophyta/Leiosporocerotopsida,Bryophytes/" CHECKME.csv

# Pteridophytes
sed -i".bak" "s/Polypodiopsida,Streptophyta/Polypodiopsida,Pteridophytes/" CHECKME.csv

# Angiosperms
sed -i".bak" "s/Magnoliopsida,Streptophyta/Magnoliopsida,Angiosperms/" CHECKME.csv
sed -i".bak" "s/Liliopsida,Streptophyta/Liliopsida,Angiosperms/" CHECKME.csv
sed -i".bak" "s/eudicots,Streptophyta/eudicots,Angiosperms/" CHECKME.csv
sed -i".bak" "s/Boraginales,NA,Streptophyta/Boraginales,NA,Angiosperms/" CHECKME.csv

```

Split Streptophyta according to http://www.theplantlist.org/

```bash
cd ~/data/organelle/plastid/summary

# Angiosperms
perl -Mojo -e '
    g(q{http://www.theplantlist.org/browse/A/})->dom
    ->find(q{li > a > i[class=family]})
    ->each( sub { print shift->text . "\n" } );
    ' > Angiosperms.tmp
echo Aceraceae >> Angiosperms.tmp
echo Asphodelaceae >> Angiosperms.tmp
echo Asteraceae >> Angiosperms.tmp
echo Campynemataceae >> Angiosperms.tmp
echo Chenopodiaceae >> Angiosperms.tmp
echo Fabaceae >> Angiosperms.tmp
echo Francoaceae >> Angiosperms.tmp
echo Hyacinthaceae >> Angiosperms.tmp
echo Nyssaceae >> Angiosperms.tmp
echo Taccaceae >> Angiosperms.tmp
echo Viscaceae >> Angiosperms.tmp

cat Angiosperms.tmp |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({},\w+,\w+),Streptophyta/\1,Angiosperms/g
        '\'' CHECKME.csv
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
        '\'' CHECKME.csv
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
        '\'' CHECKME.csv
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
3356 ---------> 3337 ---------> 2253 ---------> 2940
        NA             genus          family
```

```bash
mkdir -p ~/data/organelle/plastid/summary
cd ~/data/organelle/plastid/summary

# 3356
cat CHECKME.csv | grep -v "^#" | wc -l

# filter out accessions without linage information (strain, species, genus and family)
cat CHECKME.csv |
    perl -nla -F"," -e '
        /^#/ and next;
        ($F[2] eq q{NA} or $F[3] eq q{NA} or $F[4] eq q{NA} or $F[5] eq q{NA} ) and next;
        print
    ' \
    > valid.tmp

# 3337
wc -l valid.tmp

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

# 2253
wc -l valid.genus.tmp

#----------------------------#
# Family
#----------------------------#
# get some genera back as candidates for outgroup
cat valid.genus.tmp |
    perl -nla -F"," -e 'printf qq{,$F[5],\n}' \
    > family.tmp

# intersect between two files
grep -F -f family.tmp valid.tmp > valid.family.tmp

# 2940
wc -l valid.family.tmp

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
cd ~/data/organelle/plastid/summary

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

#Arabidopsis lyrata,2
#Astragalus mongholicus,2
#Cannabis sativa,2
#Capsicum baccatum,3
#Fragaria vesca,2
#Gossypium herbaceum,2
#Magnolia officinalis,2
#Musa balbisiana,2
#Olea europaea,4
#Oryza sativa,4
#Paris polyphylla,2
#Physcomitrella patens,2
#Pisum sativum,2
#Plasmodium falciparum,2
#Saccharum hybrid cultivar,3
#Sinalliaria limprichtiana,2
#Solanum lycopersicum,2
#Vitis aestivalis,2
#Vitis cinerea,4
#Vitis rotundifolia,2

# strain name not equal to species
cat DOWNLOAD.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e '$F[2] ne $F[3] and print $F[2]' |
    sort

#Arabidopsis lyrata subsp. lyrata
#Astragalus mongholicus var. nakaianus
#Babesia bovis T2Bo
#Babesia microti strain RI
#Brassica rapa subsp. pekinensis
#Calycanthus floridus var. glaucus
#Capsicum baccatum var. baccatum
#Capsicum baccatum var. pendulum
#Capsicum baccatum var. praetermissum
#Caragana rosea var. rosea
#Corylus ferox var. thibetica
#Cucumis melo subsp. melo
#Eucalyptus globulus subsp. globulus
#Fagopyrum esculentum subsp. ancestrale
#Fragaria vesca subsp. bracteata
#Fragaria vesca subsp. vesca
#Gossypium herbaceum subsp. africanum
#Gracilaria tenuistipitata var. liui
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
#Paris polyphylla var. chinensis
#Paris polyphylla var. yunnanensis
#Phalaenopsis aphrodite subsp. formosana
#Phyllostachys nigra var. henonis
#Pisum sativum subsp. elatius
#Plasmodium chabaudi chabaudi
#Plasmodium falciparum 3D7
#Plasmodium falciparum HB3
#Pseudotsuga sinensis var. wilsoniana
#Rosa chinensis var. spontanea
#Saccharum hybrid cultivar NCo 310
#Saccharum hybrid cultivar SP80-3280
#Sinalliaria limprichtiana var. grandifolia
#Thalassiosira oceanica CCMP1005
#Vitis aestivalis var. linsecomii
#Vitis cinerea var. cinerea
#Vitis cinerea var. floridana
#Vitis cinerea var. helleri
#Vitis rotundifolia var. munsoniana

```

Create abbreviations.

```bash
cd ~/data/organelle/plastid/summary

echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr' > ABBR.csv
cat DOWNLOAD.csv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 |
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


# count downloaded sequences # 2577 #
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

#   count every ranks
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
    --rawphylo --parallel 8 -v

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

### Batch running for groups

The old prepare_run.sh

```bash
mkdir -p ~/data/organelle/plastid.working
cd ~/data/organelle/plastid.working

bash ../plastid.cmd.txt 2>&1 | tee log_cmd.txt
# bash ../plastid.redo.cmd.txt 2>&1 | tee log_redo_cmd.txt # skip real_chr and repeatmasker

#----------------------------#
# Approach 1: one by one
#----------------------------#
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing $d <====\""
    echo bash $d/1_real_chr.sh ;
    echo bash $d/2_file_rm.sh ;
    echo bash $d/3_pair_cmd.sh ;
    echo bash $d/4_rawphylo.sh ;
    echo bash $d/5_multi_cmd.sh ;
    echo bash $d/7_multi_db_only.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

#----------------------------#
# Approach 2: step by step
#----------------------------#
# real_chr
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 1_real_chr.sh | sort `;do
    echo bash $f ;
    echo ;
done  > run_1.sh

# RepeatMasker
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 2_file_rm.sh | sort `;do
    echo bash $f ;
    echo ;
done  > run_2.sh

# pair
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 3_pair_cmd.sh | sort `;do
    echo bash $f ;
    echo ;
done  > run_3.sh

# rawphylo
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 4_rawphylo.sh | sort `;do
    echo bash $f ;
    echo ;
done  > run_4.sh

# multi cmd
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 5_multi_cmd.sh | sort `;do
    echo bash $f ;
    echo ;
done  > run_5.sh

# multi db
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 7_multi_db_only.sh | sort `;do
    echo bash $f ;
    echo ;
done  > run_7.sh

# 24 cores
cat run_1.sh | grep . | parallel -r -j 16 2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel -r -j 8  2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel -r -j 16 2>&1 | tee log_3.txt
cat run_4.sh | grep . | parallel -r -j 4  2>&1 | tee log_4.txt
cat run_5.sh | grep . | parallel -r -j 4  2>&1 | tee log_5.txt
cat run_7.sh | grep . | parallel -r -j 16 2>&1 | tee log_7.txt

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel -r rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel -r rm
```

### Self alignments.

```bash
cd ~/data/organelle/

perl -p -e '
    s/plastid\.working/plastid_self.working/g;
    s/multi_batch/self_batch/g;
    s/(\-\-parallel)/--length 1000 \1/g;
' plastid.cmd.txt > plastid_self.cmd.txt

mkdir -p ~/data/organelle/plastid_self.working
cd ~/data/organelle/plastid_self.working

time bash ../plastid_self.cmd.txt 2>&1 | tee log_cmd.txt

# Don't need 6_feature_cmd.sh
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing $d <====\""
    echo bash $d/1_real_chr.sh ;
    echo bash $d/2_file_rm.sh ;
    echo bash $d/3_self_cmd.sh ;
    echo bash $d/4_proc_cmd.sh ;
    echo bash $d/5_circos_cmd.sh ;
    echo bash $d/7_pair_stat.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

# clean
find . -mindepth 1 -maxdepth 2 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | parallel -r rm -fr

# clean mysql
#find  /usr/local/var/mysql -type d -name "[A-Z]*" | parallel -r rm -fr
```

### Alignments of families for outgroups.

```bash
mkdir -p ~/data/organelle/plastid_families
cd ~/data/organelle/plastid_families

time bash ../plastid_families.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing $d <====\""
    echo bash $d/1_real_chr.sh ;
    echo bash $d/2_file_rm.sh ;
    echo bash $d/3_pair_cmd.sh ;
    echo bash $d/4_rawphylo.sh ;
    echo bash $d/5_multi_cmd.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

find ~/data/organelle/plastid_families -type f -path "*_phylo*" -name "*.nwk"

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel -r rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel -r rm

```

In previous steps, we have manually edited `~/Scripts/withncbi/doc/plastid_OG.md` and generated
`genus_OG.tsv`.

*D* between target and outgroup should be around **0.05**.

```bash
mkdir -p ~/data/organelle/plastid_OG
cd ~/data/organelle/plastid_OG

time bash ../plastid_OG.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `; do
    echo "echo \"====> Processing $d <====\""
    echo bash $d/1_real_chr.sh ;
    echo bash $d/2_file_rm.sh ;
    echo bash $d/3_pair_cmd.sh ;
    echo bash $d/4_rawphylo.sh ;
    echo bash $d/5_multi_cmd.sh ;
    echo bash $d/7_multi_db_only.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel -r rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel -r rm

```

### LSC and SSC

IRA and IRB are presented by `plastid_self.working/${GENUS}/Results/${STRAIN}/${STRAIN}.links.tsv`.

```bash
find ~/data/organelle/plastid_self.working -type f -name "*.links.tsv" \
    | xargs wc -l | sort -n \
    | grep -v "total" \
    | perl -nl -e 's/^\s*//g; /^(\d+)\s/ and print $1' \
    | uniq -c

```

Manually check strains not containing singular link.

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
#there are 2 special strains which have no IR but a palindromic sequence.(Asa_sieboldii,Epipo_aphyllum)

Create `ir_lsc_ssc.tsv` for slicing alignments.

Manually edit it then move to `~/Scripts/withncbi/doc/ir_lsc_ssc.tsv`.

`#genus abbr    role    accession   chr_size    IR  LSC SSC`

* `NA` - not self-aligned
* `MULTI` - multiply link records
* `NONE` - no link records
* `WRONG` - unexpected link records

```bash
cd ~/data/organelle/plastid_summary

cat plastid.ABBR.csv \
    | grep -v "^#" \
    | perl -nla -F"," -MAlignDB::IntSpan -MPath::Tiny -e '
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

### Can't get clear ir information

   * Grateloupia
       *Grat_filicina
       *Grat_taiwanensis
   * Caulerpa
       *Cau_cliftonii
       *Cau_racemosa
   * Caloglossa
       *Calog_beccarii
       *Calog_intermedia
       *Calog_monosticha
   * Pisum
       *Pisum_fulvum
       *Pisum_sativum
       *Pisum_sativum_subsp_elatius
   * Dasya
       *Dasya_naccarioides
   * Diplazium
       *Diplazium_unilobum
   * Bryopsis
       *Bryop_plumosa
       *Bryop_sp_HV04063
       *Bry_sp_HV04063
   * Medicago
       *Med_falcata
       *Med_hybrida
       *Med_papillosa
       *Med_truncatula
   * Aegilops
       *Aeg_cylindrica
       *Aeg_geniculata
       *Aeg_speltoides
       *Aeg_tauschii
   * Prototheca
       *Prot_cutis
       *Prot_stagnorum
       *Prot_zopfii
   * Cryptomonas
       *Cryptomo_curvata
       *Cryptomo_paramecium
   * Monotropa
       *Monotropa_hypopitys
   * Liagora
       *Liagora_brachyclada
       *Liagora_harveyana
   * Taiwania
       *Tai_cryptomerioides
       *Tai_flousiana
   * Pinus
       *Pinus_armandii
       *Pinus_bungeana
       *Pinus_contorta
       *Pinus_gerardiana
       *Pinus_greggii
       *Pinus_jaliscana
       *Pinus_koraiensis
       *Pinus_krempfii
       *Pinus_lambertiana
       *Pinus_longaeva
       *Pinus_massoniana
       *Pinus_monophylla
       *Pinus_nelsonii
       *Pinus_oocarpa
       *Pinus_pinea
       *Pinus_sibirica
       *Pinus_strobus
       *Pinus_sylvestris
       *Pinus_tabuliformis
       *Pinus_taeda
       *Pinus_taiwanensis
       *Pinus_thunbergii
   * Taxus
       *Taxus_mairei
   * Picea
       *Pic_abies
       *Pic_asperata
       *Pic_crassifolia
       *Pic_glauca
       *Pic_jezoensis
       *Pic_morrisonicola
       *Pic_sitchensis
   * Gracilariopsis
       *Gracilario_chorda
       *Gracilario_lemaneiformis
   * Fragaria
       *Frag_mandshurica
       *Frag_vesca_subsp_bracteata
   * Ulva
       *Ulva_fasciata
       *Ulva_flexuosa
       *Ulva_linza
       *Ulva_prolifera
   * Monomorphina
       *Monom_parapyrum
   * Epipogium
       *Epipo_aphyllum
   * Euglena
       *Euglena_archaeoplastidiata
       *Euglena_viridis
   * Chlorella
       *Chlore_heliozoae
       *Chlore_sorokiniana
       *Chlore_variabilis
       *Chlore_vulgaris
   * Erodium
       *Ero_carvifolium
       *Ero_crassifolium
       *Ero_manescavi
       *Ero_rupestre
   * Larix
       *Lar_decidua
       *Lar_sibirica
   * Amentotaxus
       *Ame_argotaenia
       *Ame_formosana
   * Pyropia
       *Pyro_perforata
   * Ceramium
       *Ceram_cimbricum
       *Ceram_japonicum
       *Ceram_sungminbooi
   * Hildenbrandia
       *Hilde_rivularis
       *Hilde_rubra
   * Pilostyles
       *Pilo_aethiopica
       *Pilo_hamiltonii
   * Codium
       *Codi_decorticatum
       *Codi_sp_arenicola
   * Torreya
       *Torreya_fargesii
       *Torreya_grandis
   * Vertebrata
       *Vert_australis
       *Vert_isogona
       *Vert_lanosa
       *Vert_thuyoides
   * Wisteria
       *Wis_floribunda
       *Wis_sinensis
   * Phelipanche
       *Pheli_purpurea
       *Pheli_ramosa
   * Glycyrrhiza
       *Glycy_glabra
       *Glycy_lepidota
   * Cephalotaxus
       *Cephalo_wilsoniana
   * Polysiphonia
       *Polysi_brodiei
       *Polysi_elongata
       *Polysi_infestans
       *Polysi_schneideri
       *Polysi_scopulorum
       *Polysi_sertularioides
       *Polysi_stricta
   * Lathyrus
       *Lathy_clymenum
       *Lathy_davidii
       *Lathy_graminifolius
       *Lathy_inconspicuus
       *Lathy_littoralis
       *Lathy_ochroleucus
       *Lathy_odoratus
       *Lathy_palustris
       *Lathy_sativus
       *Lathy_tingitanus
   * Babesia
       *Bab_orientalis
   * Gelidium
       *Gelidi_elegans
       *Gelidi_vagum
   * Bostrychia
       *Bos_moritziana
       *Bos_simpliciuscula
       *Bos_tenella
   * Membranoptera
       *Mem_platyphylla
       *Mem_tenuis
       *Mem_weeksiae
   * Astragalus
       *Astra_mongholicus
       *Astra_mongholicus_var_nakaianus
   * Asarum
       *Asa_minus
       *Asa_sieboldii
   * Gracilaria
       *Gracilaria_changii
       *Gracilaria_chilensis
       *Gracilaria_firma
       *Gracilaria_salicornia
       *Gracilaria_tenuistipitata_var_liui
       *Gracilaria_vermiculophylla
   * Plasmodium
       *Plas_chabaudi_chabaudi
       *Plas_falciparum_HB3
       *Plas_gallinaceum
       *Plas_relictum
       *Plas_vivax
   * Triticum
       *Trit_urartu
   * Trifolium
       *Trif_aureum
       *Trif_boissieri
       *Trif_glanduliferum
       *Trif_grandiflorum
       *Trif_strictum

#*:there are 2 special strains which has only one palindromic sequence rather than IR.(as mentioned before)

### Slices of IR, LSC and SSC

Without outgroups.

Be cautious to alignments with low coverage.

```bash
mkdir -p ~/data/organelle/plastid_slices
cd ~/data/organelle/plastid_slices

cat ~/Scripts/withncbi/doc/ir_lsc_ssc.tsv \
    | perl -nla -F"\t" -MAlignDB::IntSpan -Mstrict -Mwarnings -e '
        /^#/ and next;
        $F[2] eq q{Target} or next;
        $F[5] =~ /\d+/ or next;

        print qq{# $F[0]};

        next unless -e "$ENV{HOME}/data/organelle/plastid.working/$F[0]/$F[0]_refined/$F[3].synNet.maf.gz.fas.gz";

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
            print qq{    ~/data/organelle/plastid.working/$F[0]/$F[0]_refined/$F[3].synNet.maf.gz.fas.gz \\};
            print qq{    $F[0].$key.yml};
            print qq{perl ~/Scripts/alignDB/alignDB.pl \\};
            print qq{    -d $F[0]_$key \\};
            print qq{    -da ~/data/organelle/plastid_slices/$F[0].$key.fas \\};
            print qq{    -a ~/data/organelle/plastid.working/$F[0]/Stats/anno.yml \\};
            print qq{    -chr ~/data/organelle/plastid.working/$F[0]/chr_length.csv \\};
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

#### SNP t-test

```bash
mkdir -p ~/data/organelle/plastid_summary/ttest
cd ~/data/organelle/plastid_summary/ttest

cat ~/Scripts/withncbi/doc/plastid_OG.md \
    | grep -v "^#" | grep . \
    | cut -d',' -f 1 \
    | parallel echo perl ~/Scripts/alignDB/stat/multi_stat_factory.pl -d {}_OG -r 52 \
    > plastid_og_ttest.cmd.txt

sh plastid_og_ttest.cmd.txt
```

## Cyanobacteria

### Genus and Species counts

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

### Use `bac_prepare.pl`

```bash
mkdir -p ~/data/organelle/cyanobacteria
cd ~/data/organelle/cyanobacteria

perl ~/Scripts/withncbi/taxon/bac_prepare.pl --db gr_prok \
    --seq_dir ~/data/organelle/cyanobacteria/bac_seq_dir \
    -p 1218 --get_seq -n Prochlorococcus

perl ~/Scripts/withncbi/taxon/bac_prepare.pl --db gr_prok \
    --seq_dir ~/data/organelle/cyanobacteria/bac_seq_dir \
    -p 74546,93060,146891,167546 --get_seq -t 74546 -n Prochlorococcus_marinus

perl ~/Scripts/withncbi/taxon/bac_prepare.pl --db gr_prok \
    --seq_dir ~/data/organelle/cyanobacteria/bac_seq_dir \
    -p 74546,93060,146891,167546,59919 --get_seq -t 74546 -o 59919 -n Prochlorococcus_marinus_OG

perl ~/Scripts/withncbi/taxon/bac_prepare.pl --db gr_prok \
    --seq_dir ~/data/organelle/cyanobacteria/bac_seq_dir \
    -p 1142 --get_seq -n Synechocystis

perl ~/Scripts/withncbi/taxon/bac_prepare.pl --db gr_prok \
    --seq_dir ~/data/organelle/cyanobacteria/bac_seq_dir \
    -p 1148 --get_seq -n Synechocystis_sp_PCC_6803

perl ~/Scripts/withncbi/taxon/bac_prepare.pl --db gr_prok \
    --seq_dir ~/data/organelle/cyanobacteria/bac_seq_dir \
    -p 1148,1147 -o 1147 --get_seq -n Synechocystis_sp_PCC_6803_OG

export BAC_DIR=Prochlorococcus
sh ${BAC_DIR}/prepare.sh
sh ${BAC_DIR}/1_real_chr.sh
sh ${BAC_DIR}/2_file_rm.sh
sh ${BAC_DIR}/3_pair_cmd.sh
sh ${BAC_DIR}/4_rawphylo.sh
sh ${BAC_DIR}/5_multi_cmd.sh
sh ${BAC_DIR}/7_multi_db_only.sh
unset BAC_DIR

sh Prochlorococcus/prepare.sh
sh Prochlorococcus_marinus/prepare.sh
sh Prochlorococcus_marinus_OG/prepare.sh
sh Synechocystis/prepare.sh
sh Synechocystis_sp_PCC_6803/prepare.sh
sh Synechocystis_sp_PCC_6803_OG/prepare.sh

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d -not -path "*bac_seq_dir" | sort `;do \
    echo bash $d/1_real_chr.sh ; \
    echo bash $d/2_file_rm.sh ; \
    echo bash $d/3_pair_cmd.sh ; \
    echo bash $d/4_rawphylo.sh ; \
    echo bash $d/5_multi_cmd.sh ; \
    echo bash $d/7_multi_db_only.sh ; \
    echo ; \
done  > runall.sh


find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel -r rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel -r rm
```

## Summary

### Copy xlsx files

```bash
mkdir -p ~/data/organelle/plastid_summary/xlsx
cd ~/data/organelle/plastid_summary/xlsx

find  ~/data/organelle/plastid.working -type f -name "*.common.xlsx" \
    | grep -v "vs[A-Z]" \
    | parallel cp {} .

find  ~/data/organelle/plastid_OG -type f -name "*.common.xlsx" \
    | grep -v "vs[A-Z]" \
    | parallel cp {} .

```

### Genome list

Create `plastid.list.csv` from `plastid.GENUS.csv` with sequence
lengths.

```bash
mkdir -p ~/data/organelle/plastid_summary/table
cd ~/data/organelle/plastid_summary/table

# manually set orders in `plastid_OG.md`
echo "#genus" > genus_all.lst
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
    >> genus_all.lst

echo "#abbr,genus,accession,length" > length.tmp
find ~/data/organelle/plastid.working -type f -name "chr.sizes" \
    | parallel --jobs 1 --keep-order -r '
        perl -nl -e '\''
            BEGIN {
                %l = ();
            }
            
            next unless /\w+\t\d+/;
            my ($key, $value) = split /\t/;
            $l{$key} = $value;
            
            END {
                my $chrs = join "|", sort keys %l;
                my $length = 0;
                $length += $_ for values %l;
                
                $ARGV =~ /working\/(\w+)\/(\w+)\/(\w+)/;
                print qq{$3,$1,$chrs,$length}
            }
        '\'' \
        {}
    ' \
    >> length.tmp

echo "#abbr,phylum,family,genus,taxon_id" > abbr.tmp
cat ~/data/organelle/plastid_summary/plastid.GENUS.csv \
    | grep -v "^#" \
    | perl -nla -F"," -e 'print qq{$F[9],$F[8],$F[5],$F[4],$F[0]}' \
    >> abbr.tmp

# #abbr,genus,accession,length,phylum,family,genus,taxon_id
cat length.tmp abbr.tmp \
    | perl ~/Scripts/withncbi/util/merge_csv.pl \
        -f 0 --concat -o stdout \
    | perl -nl -a -F"," -e 'print qq{$F[4],$F[5],$F[6],$F[0],$F[7],$F[2],$F[3]}' \
    > list.tmp

echo "#phylum,family,genus,abbr,taxon_id,accession,length" > plastid.list.csv
cat list.tmp \
    | grep -v "#" \
    | perl -nl -a -F',' -MPath::Tiny -e '
        BEGIN{
            %genus, %target;
            my @l1 = path(q{genus_all.lst})->lines({ chomp => 1});
            $genus{$l1[$_]} = $_ for (0 .. $#l1);
        }
        my $idx = $genus{$F[2]};
        die qq{$_\n} unless defined $idx;
        print qq{$_,$idx};
    ' \
    | sort -n -t',' -k8,8 \
    | cut -d',' -f 1-7 \
    >> plastid.list.csv
# 1905

rm *.tmp
```

### Genome alignment statistics

Some genera will be filtered out here.

Criteria:

* Coverage >= 0.4
* Total number of indels >= 100
* Genome D < 0.2

```bash
mkdir -p ~/data/organelle/plastid_summary/table

cd ~/data/organelle/plastid_summary/xlsx
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

cat ~/data/organelle/plastid_summary/table/genus_all.lst \
    | grep -v "^#" \
    | TT_FILE=Table_alignment.tt perl -MTemplate -nl -e '
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

cp -f Table_alignment_all.xlsx ~/data/organelle/plastid_summary/table
cp -f Table_alignment_all.csv ~/data/organelle/plastid_summary/table

cd ~/data/organelle/plastid_summary/table

echo "Genus,avg_size" > group_avg_size.csv
cat plastid.list.csv \
    | grep -v "#" \
    | perl -nla -F"," -e '
        $count{$F[2]}++;
        $sum{$F[2]} += $F[6];
        END {
            for $k (sort keys %count) {
                printf qq{%s,%d\n}, $k, $sum{$k}/$count{$k};
            }
        }
    ' \
    >> group_avg_size.csv

cat Table_alignment_all.csv group_avg_size.csv \
    | perl ~/Scripts/withncbi/util/merge_csv.pl \
    -f 0 --concat -o stdout \
    > Table_alignment_all.1.csv

echo "Genus,coverage" > group_coverage.csv
cat Table_alignment_all.1.csv \
    | perl -nla -F',' -e '
        $F[2] =~ /[\.\d]+/ or next;
        $F[6] =~ /[\.\d]+/ or next;
        $c = $F[2] * 1000 * 1000 / $F[6];
        print qq{$F[0],$c};
    ' \
    >> group_coverage.csv

cat Table_alignment_all.1.csv group_coverage.csv \
    | perl ~/Scripts/withncbi/util/merge_csv.pl \
    -f 0 --concat -o stdout \
    > Table_alignment_all.2.csv

echo "Genus,indels" > group_indels.csv
cat Table_alignment_all.2.csv \
    | perl -nla -F',' -e '
        $F[6] =~ /[\.\d]+/ or next;
        $c = $F[3] / 100 * $F[2] * 1000 * 1000;
        print qq{$F[0],$c};
    ' \
    >> group_indels.csv

cat Table_alignment_all.2.csv group_indels.csv \
    | perl ~/Scripts/withncbi/util/merge_csv.pl \
    -f 0 --concat -o stdout \
    > Table_alignment_for_filter.csv

# real filter
cat Table_alignment_for_filter.csv \
    | perl -nla -F',' -e '
        $F[6] =~ /[\.\d]+/ or next;
        $F[0] =~ s/"//g;
        print $F[0] if ($F[7] < 0.4 or $F[8] < 100 or $F[4] > 0.2);
    ' \
    > genus_exclude.lst

grep -v -Fx -f genus_exclude.lst genus_all.lst > genus.lst

rm ~/data/organelle/plastid_summary/table/Table_alignment_all.[0-9].csv
rm ~/data/organelle/plastid_summary/table/group_*csv

#
cd ~/data/organelle/plastid_summary/xlsx
cat ~/data/organelle/plastid_summary/table/genus.lst \
    | grep -v "^#" \
    | TT_FILE=Table_alignment.tt perl -MTemplate -nl -e '
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

cp -f ~/data/organelle/plastid_summary/xlsx/Table_alignment.xlsx ~/data/organelle/plastid_summary/table
cp -f ~/data/organelle/plastid_summary/xlsx/Table_alignment.csv ~/data/organelle/plastid_summary/table
```

### Groups

```bash
mkdir -p ~/data/organelle/plastid_summary/group
cd ~/data/organelle/plastid_summary/group

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

grep -Fx -f ~/data/organelle/plastid_summary/table/genus.lst Angiosperms.txt > Angiosperms.lst
grep -Fx -f ~/data/organelle/plastid_summary/table/genus.lst Gymnosperms.txt > Gymnosperms.lst

find . -type f -name "*.txt" \
    | xargs cat \
    | grep -Fx -f ~/data/organelle/plastid_summary/table/genus.lst \
    > Others.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nla -F',' -e '$F[1] > 0.05 and print $F[0];' \
    > group_4.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nla -F',' -e '$F[1] > 0.02 and $F[1] <= 0.05 and print $F[0];' \
    > group_3.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nla -F',' -e '$F[1] > 0.005 and $F[1] <= 0.02 and print $F[0];' \
    > group_2.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nla -F',' -e '$F[1] <= 0.005 and print $F[0];' \
    > group_1.lst

rm *.txt
```

NCBI Taxonomy tree

```bash
cd ~/data/organelle/plastid_summary/group

cat ~/data/organelle/plastid_summary/table/genus.lst \
    | grep -v "^#" \
    | perl -e '
        @ls = <>;
        $str = qq{bp_taxonomy2tree.pl \\\n};
        for (@ls) {
            chomp;
            $str .= qq{    -s "$_" \\\n};
        }
        $str .= qq{    -e \n};
        print $str
    ' \
    > genera_tree.sh

bash genera_tree.sh > genera.tree
```

### Phylogenic trees of each genus with outgroup

```bash
mkdir -p ~/data/organelle/plastid_summary/trees

cat ~/Scripts/withncbi/doc/plastid_OG.md \
    | grep -v "^#" \
    | grep . \
    | cut -d',' -f 1 \
    > ~/data/organelle/plastid_summary/trees/list.txt

find ~/data/organelle/plastid_OG -type f -path "*_phylo*" -name "*.nwk" \
    | parallel -j 1 cp {} ~/data/organelle/plastid_summary/trees
```

### d1, d2

`collect_xlsx.pl`

```bash
cd ~/data/organelle/plastid_summary/xlsx

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

cat ~/data/organelle/plastid_summary/table/genus.lst \
    | grep -v "^#" \
    | TT_FILE=cmd_collect_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, }) 
                or die Template->error;
        }
    ' \
    > cmd_collect_d1_d2.sh

cd ~/data/organelle/plastid_summary/xlsx
bash cmd_collect_d1_d2.sh
```

`sep_chart.pl`

```bash
mkdir -p ~/data/organelle/plastid_summary/fig

cd ~/data/organelle/plastid_summary/xlsx

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

cat ~/data/organelle/plastid_summary/group/group_1.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
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

cat ~/data/organelle/plastid_summary/group/group_2.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.03,
                y_max2 => 0.04,
                postfix => q{group_2}, }) 
                or die Template->error;
        }
    ' \
    > cmd_chart_group_2.sh

cat ~/data/organelle/plastid_summary/group/group_3.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
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

cat ~/data/organelle/plastid_summary/group/group_4.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.15,
                y_max2 => 0.15,
                postfix => q{group_4}, }) 
                or die Template->error;
        }
    ' \
    > cmd_chart_group_4.sh

bash cmd_chart_group_1.sh
bash cmd_chart_group_2.sh
bash cmd_chart_group_3.sh
bash cmd_chart_group_4.sh

rm ~/data/organelle/plastid_summary/xlsx/*.csv
cp ~/data/organelle/plastid_summary/xlsx/*.pdf ~/data/organelle/plastid_summary/fig

# Coreldraw doesn't play well with computer modern fonts (latex math).
# perl ~/Scripts/fig_table/tikz_chart.pl -i cmd_plastid_d1_A2A8_B2B8.group_1.csv -xl 'Distance to indels ($d_1$)' -yl 'Nucleotide divergence ($D$)' --y_min 0.0 --y_max 0.01 -x_min 0 -x_max 5 --style_dot --pdf
```
