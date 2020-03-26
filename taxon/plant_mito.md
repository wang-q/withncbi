# Process plant mitochondrion genomes

[TOC levels=1-3]: # ""

- [Process plant mitochondrion genomes](#process-plant-mitochondrion-genomes)
- [Update taxdmp](#update-taxdmp)
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
  - [Restrict taxonomy ids to green plants with `taxon/id_restrict.pl`.](#restrict-taxonomy-ids-to-green-plants-with-taxonid_restrictpl)
- [Add lineage information](#add-lineage-information)
  - [Can't get clear taxonomy](#cant-get-clear-taxonomy)
- [Filtering based on valid families and genera](#filtering-based-on-valid-families-and-genera)
- [Find a way to name these](#find-a-way-to-name-these)
- [Download sequences and regenerate lineage information](#download-sequences-and-regenerate-lineage-information)
  - [Numbers for higher ranks](#numbers-for-higher-ranks)
  - [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
- [Prepare sequences for lastz](#prepare-sequences-for-lastz)
- [Aligning without outgroups](#aligning-without-outgroups)
  - [Create `mito_t_o.md` for picking targets and outgroups.](#create-mito_t_omd-for-picking-targets-and-outgroups)
  - [Create alignments plans without outgroups](#create-alignments-plans-without-outgroups)
  - [Plans for align-able targets](#plans-for-align-able-targets)
  - [Aligning w/o outgroups](#aligning-wo-outgroups)
  - [Aligning with outgroups](#aligning-with-outgroups)
- [Self alignments](#self-alignments)
- [Summary](#summary)
  - [Copy xlsx files](#copy-xlsx-files)
  - [Genome list](#genome-list)
  - [Statistics of genome alignments](#statistics-of-genome-alignments)
  - [Groups](#groups)
  - [Phylogenic trees of each genus with outgroup](#phylogenic-trees-of-each-genus-with-outgroup)
  - [d1, d2](#d1-d2)


# Update taxdmp

*Update `~/data/NCBI/taxdmp` before running `id_restrict.pl` or `id_project_to.pl`*.

# Scrap id and acc from NCBI

Open browser and visit
[NCBI mitochondrion page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=organelle).
Save page to a local file, html only. In this case, it's
`doc/green_plants_mitochondrion_181207.html`.

All [Eukaryota](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=organelle),
`doc/eukaryota_mitochondrion_181207.html`.

```text
Eukaryota (2759)                8746
    Viridiplantae (33090)       221
        Chlorophyta (3041)      51
        Streptophyta (35493)    170
```

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions from all history
pages.

Got **10277** accessions.

```bash
mkdir -p ~/data/organelle/mito/GENOMES
cd ~/data/organelle/mito/GENOMES

rm webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/eukaryota_mitochondrion_181207.html \
    >> webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/green_plants_mitochondrion_181207.html \
    >> webpage_id_seq.csv

```

Use `taxon/gb_taxon_locus.pl` to extract information from refseq genbank files.

```bash
cd ~/data/organelle/mito/GENOMES

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz

gzip -dcf mitochondrion.*.genomic.gbff.gz > genomic.gbff

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv

rm genomic.gbff

cat refseq_id_seq.csv | grep -v "^#" | wc -l
# 10277

# combine
cat webpage_id_seq.csv refseq_id_seq.csv |
    sort -u | # duplicated id-seq pair
    sort -t, -k1,1 \
    > id_seq.csv

cat id_seq.csv | grep -v "^#" | wc -l
# 10277

```

## Restrict taxonomy ids to green plants with `taxon/id_restrict.pl`.

```bash
cd ~/data/organelle/mito/GENOMES

# Viridiplantae 33090
echo '#strain_taxon_id,accession' > plant_id_seq.csv
cat id_seq.csv |
    grep -v "^#" |
    perl ~/Scripts/withncbi/taxon/id_restrict.pl -s "," -a 33090 \
    >> plant_id_seq.csv

# 293
cat plant_id_seq.csv | grep -v "^#" | wc -l

cat plant_id_seq.csv |
    cut -d',' -f 1 |
    sort -n |
    uniq -c |
    grep -v -E '\s+1\s+'
#      3 3659 Cucumis sativus has 3 chromosomes
#      2 3708 Brassica napus linear plasmid NC_004946
#      2 39946 Oryza sativa indica plasmid B2 NC_001776
#      2 39947 Oryza sativa japonica plasmid B1 NC_001751
#      2 51329 Polytomella parva has 2 chromosomes
#      2 351366 Polytomella sp. SAG 63-10 has 2 chromosomes

#      2 3702 Arabidopsis thaliana NC_001284.2 was removed by RefSeq staff

sed -i".bak" "/,NC_001284$/d" plant_id_seq.csv # Arabidopsis thaliana
sed -i".bak" "/,NC_004946$/d" plant_id_seq.csv # Brassica napus
sed -i".bak" "/,NC_001751$/d" plant_id_seq.csv # Oryza sativa japonica
sed -i".bak" "/,NC_001776$/d" plant_id_seq.csv # Oryza sativa indica

# Vicia faba mitochondrial plasmid MtVFPL3 NC_011084
# Zea mays mitochondrial plasmid pBMSmt1.9 NC_001400
sed -i".bak" "/,NC_011084$/d" plant_id_seq.csv
sed -i".bak" "/,NC_001400$/d" plant_id_seq.csv

```

# Add lineage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `CHECKME.csv`.

```bash
mkdir -p ~/data/organelle/mito/summary
cd ~/data/organelle/mito/summary

# generate a .csv file for manually checking
echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum' > CHECKME.csv
cat ../GENOMES/plant_id_seq.csv |
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

```bash
cd ~/data/organelle/mito/summary

# darwin (bsd) need "" for -i
sed -i".bak" "s/\'//g" CHECKME.csv

# Anomodon attenuatus and Anomodon rugelii (Bryophytes) was grouped to Streptophyta.
sed -i".bak" "s/Bryopsida,Streptophyta/Bryopsida,Bryophytes/" CHECKME.csv

sed -i".bak" "s/Pycnococcaceae,NA/Pycnococcaceae,Pseudoscourfieldiales/" CHECKME.csv

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
cd ~/data/organelle/mito/summary

#perl -Mojo -e '
#    g(q{http://www.theplantlist.org/browse/A/})->dom
#    ->find(q{li > a > i[class=family]})
#    ->each( sub { print shift->text . "\n" } );
#    ' > Angiosperms.tmp

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

## Can't get clear taxonomy

* Species
  + Trebouxiophyceae sp. MX-AZ01


# Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

```text
289 ---------> 286 ---------> 113 ---------> 158
        NA           genus         family
```


```bash
mkdir -p ~/data/organelle/mito/summary
cd ~/data/organelle/mito/summary

cat CHECKME.csv | grep -v "^#" | wc -l
# 289

# filter out accessions without linage information (strain, species, genus and family)
cat CHECKME.csv |
    perl -nla -F"," -e '
        /^#/ and next;
        ($F[2] eq q{NA} or $F[3] eq q{NA} or $F[4] eq q{NA} or $F[5] eq q{NA} ) and next;
        print
    ' \
    > valid.tmp

wc -l valid.tmp
# 286

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
# 113

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
# 158

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
cd ~/data/organelle/mito/summary

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
#Beta vulgaris,2
#Cucumis sativus,3
#Oryza sativa,2
#Polytomella parva,2
#Polytomella piriformis,2
#Zea mays,2

# strain name not equal to species
cat DOWNLOAD.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e '$F[2] ne $F[3] and print $F[2]' |
    sort
#Aegilops speltoides var. ligustica
#Beta vulgaris subsp. maritima
#Beta vulgaris subsp. vulgaris
#Brassica rapa subsp. oleifera
#Marchantia polymorpha subsp. ruderalis
#Oryza sativa Indica Group
#Oryza sativa Japonica Group
#Zea mays subsp. mays
#Zea mays subsp. parviglumis

```

Create abbreviations.

```bash
cd ~/data/organelle/mito/summary

echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr' > ABBR.csv
cat DOWNLOAD.csv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 --shortsub |
    sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> ABBR.csv

```

# Download sequences and regenerate lineage information

```bash
cd ~/data/organelle/mito/GENOMES

echo "#strain_name,accession,strain_taxon_id" > name_acc_id.csv
cat ../summary/ABBR.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e 'print qq{$F[9],$F[1],$F[0]}' |
    sort \
    >> name_acc_id.csv

# Local, Runtime 4 seconds.
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

## Numbers for higher ranks

21 orders, 22 families, 36 genera and 106 species.

```bash
cd ~/data/organelle/mito/summary

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

wc -l GENUS.csv
# 113

# count every ranks
cut -d',' -f 4 GENUS.csv | sort | uniq > species.list.tmp
cut -d',' -f 5 GENUS.csv | sort | uniq > genus.list.tmp
cut -d',' -f 6 GENUS.csv | sort | uniq > family.list.tmp
cut -d',' -f 7 GENUS.csv | sort | uniq > order.list.tmp
wc -l order.list.tmp family.list.tmp genus.list.tmp species.list.tmp
#  21 order.list.tmp
#  22 family.list.tmp
#  36 genus.list.tmp
# 106 species.list.tmp

# create again with headers
grep -F -f genus.tmp ABBR.csv > GENUS.tmp

# sort by multiply columns, phylum, order, family, genus, accession
head -n 1 ABBR.csv > GENUS.csv
cat GENUS.tmp |
    sort -t',' -k9,9 -k7,7 -k6,6 -k5,5 -k2,2 \
    >> GENUS.csv

# clean
rm *.tmp *.bak

```

## Raw phylogenetic tree by MinHash

```bash
mkdir -p ~/data/organelle/mito/mash
cd ~/data/organelle/mito/mash

for name in $(cat ../summary/ABBR.csv | sed -e '1d' | cut -d"," -f 10 | sort); do
    2>&1 echo "==> ${name}"

    if [[ -e ${name}.msh ]]; then
        continue
    fi

    find ../GENOMES/${name} -name "*.fa" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${name}" -o ${name}
done

cd ~/data/organelle/mito/summary
mash triangle -E -p 8 -l <(
        cat ABBR.csv |
            sed '1d' |
            cut -d"," -f 10 |
            parallel echo "../mash/{}.msh"
    ) \
    > dist.tsv

# fill matrix with lower triangle
tsv-select -f 1-3 dist.tsv |
    (tsv-select -f 2,1,3 dist.tsv && cat) |
    (
        cut -f 1 dist.tsv |
            tsv-uniq |
            parallel -j 1 --keep-order 'echo -e "{}\t{}\t0"' &&
        cat
    ) \
    > dist_full.tsv

cat dist_full.tsv |
    Rscript -e '
        library(readr);
        library(tidyr);
        library(ape);
        pair_dist <- read_tsv(file("stdin"), col_names=F);
        tmp <- pair_dist %>%
            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0), values_fn = list(X3 = mean) )
        tmp <- as.matrix(tmp)
        mat <- tmp[,-1]
        rownames(mat) <- tmp[,1]

        dist_mat <- as.dist(mat)
        clusters <- hclust(dist_mat, method = "ward.D2")
        tree <- as.phylo(clusters)
        write.tree(phy=tree, file="tree.nwk")

        group <- cutree(clusters, h=0.3) # k=3
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

nw_display -s -b 'visibility:hidden' -w 600 -v 30 tree.nwk |
    rsvg-convert -o plant_mito.png

```

# Prepare sequences for lastz

```bash
cd ~/data/organelle/mito/GENOMES

find . -maxdepth 1 -mindepth 1 -type d |
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

## Create `mito_t_o.md` for picking targets and outgroups.

Manually edit it then move to `~/Scripts/withncbi/doc/mito_t_o.md`.

```bash
cd ~/data/organelle/mito/summary

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
            printf qq{%s,%s\n}, $genus, $F[9];
        }
    ' \
    > mito_t_o.md

# 36
cat mito_t_o.md |
    grep -v '^#' |
    grep -E '\S+' |
    wc -l

# mv mito_t_o.md ~/Scripts/withncbi/doc/mito_t_o.md

```

## Create alignments plans without outgroups

```text
GENUS.csv
#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr
```

**36** genera and **15** families.

```bash
mkdir -p ~/data/organelle/mito/taxon
cd ~/data/organelle/mito/taxon

echo -e "#Serial\tGroup\tCount\tTarget" > group_target.tsv

cat ../summary/GENUS.csv |
    grep -v "^#" |
    SERIAL=1 perl -na -F"," -MPath::Tiny -e '
        BEGIN{
            $name = q{};
            %id_of = ();
            %h = ();
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{~/Scripts/withncbi/doc/mito_t_o.md})->lines({chomp => 1});
            for (@ls) {
                @fs = split(/,/);
                $h{$fs[0]}= $fs[1];
            }
            undef @ls;
        }

        chomp for @F;
        $F[4] =~ s/\W+/_/g;
        if ($F[4] ne $name) {
            if ($name) {
                if (exists $h{$name}) {
                    my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
                    my $t = $h{$name};
                    printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
                    path(qq{$name})->spew(map {qq{$_\n}} @s);
                    $ENV{SERIAL}++;
                }
            }
            $name = $F[4];
            %id_of = ();
        }
        $id_of{$F[9]} = $F[0]; # same strain multiple chromosomes collapsed here

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = $h{$name};
            printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
            path(qq{$name})->spew(map {qq{$_\n}} @s);
        }' \
    >> group_target.tsv

cat ../summary/GENUS.csv |
    grep -v "^#" |
    SERIAL=501 perl -na -F"," -MPath::Tiny -e '
        BEGIN{
            our $name = q{};
            our %id_of = ();
        }

        chomp for @F;
        my $family = $F[5];
        $family =~ s/\W+/_/g;
        if ($family ne $name) {
            if ($name) {
                # sort by taxonomy_id
                my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
                my $t = $s[0];
                if (scalar @s > 2) {
                    printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
                    path(qq{$name})->spew(map {qq{$_\n}} @s);
                    $ENV{SERIAL}++;
                }
            }
            $name = $family;
            %id_of = ();
        }
        $id_of{$F[9]} = $F[0]; # multiple chromosomes collapsed here

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = $s[0];
            if (scalar @s > 2) {
                printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
                path(qq{$name})->spew(map {qq{$_\n}} @s);
            }
        }
    '  \
    >> group_target.tsv

```

## Plans for align-able targets

```bash
cd ~/data/organelle/mito/taxon

cat ~/Scripts/withncbi/doc/mito_t_o.md |
    grep -v "^#" |
    grep . |
    perl -nla -F"," -e 'print $F[1] ' \
    > targets.tmp

cat ../summary/groups.tsv |
    grep -F -w -f targets.tmp |
    perl -nla -F"\t" -e '($g, $s) = split q{_}, $F[1]; print qq{$F[0]\t${g}_${s}}' |
    tsv-summarize --group-by 1 --count |
    tsv-filter --ge 2:2 |
    cut -f 1 \
    > groups.tmp

cat ../summary/groups.tsv |
    grep -F -w -f targets.tmp |
    grep -F -w -f groups.tmp |
    tsv-summarize --group-by 1 --values 2 \
    > groups.lst.tmp

cat groups.lst.tmp |
    grep -v "^#" |
    SERIAL=901 perl -na -F"\t" -MPath::Tiny -e '
        chomp for @F;
        my $group = $F[0];
        $group = "group_${group}";
        my @targets = split /\|/, $F[1];

        printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $group, scalar @targets, $targets[0];
        path(qq{$group})->spew(map {qq{$_\n}} @targets);
        $ENV{SERIAL}++;
    '  \
    >> group_target.tsv

rm *.tmp

```

## Aligning w/o outgroups

```bash
cd ~/data/organelle/mito/

# genus
cat taxon/group_target.tsv |
    tsv-filter -H  --ge 1:1 --le 1:500 |
    sed -e '1d' | #grep "^200" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        if bjobs -w | grep -w {2}; then
            exit;
        fi

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/genus/{2} \
            --rawphylo --parallel 4 -v

        bash groups/genus/{2}/1_pair.sh
        bash groups/genus/{2}/2_rawphylo.sh
        bash groups/genus/{2}/3_multi.sh
    '

# family
cat taxon/group_target.tsv |
    tsv-filter -H --ge 1:501 --le 1:900 |
    sed -e '1d' | #grep "^503" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/family/{2} \
            --rawphylo --parallel 4 -v

        bash groups/family/{2}/1_pair.sh
        bash groups/family/{2}/2_rawphylo.sh
        bash groups/family/{2}/3_multi.sh
    '

# group
cat taxon/group_target.tsv |
    tsv-filter -H --ge 1:901 |
    sed -e '1d' | #grep "^915" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/group/{2} \
            --rawphylo --parallel 4 -v

        bash groups/group/{2}/1_pair.sh
        bash groups/group/{2}/2_rawphylo.sh
        bash groups/group/{2}/3_multi.sh
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

# check status
echo \
    $(find groups/genus -mindepth 1 -maxdepth 1 -type d | wc -l) \
    $(find groups/genus -mindepth 1 -maxdepth 3 -type f -name "*.nwk.pdf" | grep -w raw -v | wc -l)

find groups/genus -mindepth 1 -maxdepth 1 -type d |
    parallel -j 4 '
        lines=$(find {} -type f -name "*.nwk.pdf" | grep -w raw -v | wc -l)
        if [ $lines -eq 0 ]; then
            lines=$(find {} -type d -name "mafSynNet" | wc -l)
            if [ $lines -gt 2 ]; then
                echo {}
            fi
        fi
    ' |
    sort

```

## Aligning with outgroups

* Review alignments and phylogenetic trees generated in `groups/family/` and `groups/group/`

* Add outgroups to `mito_t_o.md` manually.

* *D* between target and outgroup should be around **0.05**.

```bash
cd ~/data/organelle/mito/

# species_og
cat taxon/group_target.tsv |
    tsv-filter -H --le 1:500 |
    sed -e '1d' | #grep "^163" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        outgroup=$(
            cat ~/Scripts/withncbi/doc/mito_t_o.md |
                grep -v "^#" |
                grep . |
                perl -nl -e '\'' m/{2},{4},(\w+)/ and print $1 '\''
            )

        if [ "${outgroup}" = "" ]; then
            exit;
        fi

        if [ ! -d "GENOMES/${outgroup}" ]; then
            exit;
        fi

        echo -e "==> Group: [{2}]\tTarget: [{4}]\tOutgroup: [${outgroup}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            GENOMES/${outgroup} \
            --multi -o groups/genus_og/{2} \
            --outgroup ${outgroup} \
            --rawphylo --parallel 4 -v

        bash groups/genus_og/{2}/1_pair.sh
        bash groups/genus_og/{2}/2_rawphylo.sh
        bash groups/genus_og/{2}/3_multi.sh
    '

```

# Self alignments

```bash
cd ~/data/organelle/mito/summary

cat <<'EOF' > egaz_templates_self.tt

# [% name %]
egaz template \
    ~/data/organelle/mito/GENOMES/[% t %] \
[% FOREACH q IN qs -%]
    ~/data/organelle/mito/GENOMES/[% q %] \
[% END -%]
    --self -o [% name %] \
    --taxon ~/data/organelle/mito/GENOMES/taxon_ncbi.csv \
    --circos --parallel 8 -v

EOF

# every genera
echo "mkdir -p ~/data/organelle/mito/self"  > ../self.cmd.txt
echo "cd       ~/data/organelle/mito/self" >> ../self.cmd.txt
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
mkdir -p ~/data/organelle/mito/self
cd ~/data/organelle/mito/self

time bash ../self.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing ${d} <====\""
    echo bash ${d}/1_self.sh;
    echo bash ${d}/3_proc.sh;
    echo bash ${d}/4_circos.sh;
    echo bash ${d}/6_chr_length.sh;
    echo bash ${d}/7_self_aligndb.sh;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

```

# Summary

## Copy xlsx files

```bash
mkdir -p ~/data/organelle/mito/summary/xlsx
cd ~/data/organelle/mito/summary/xlsx

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
mkdir -p ~/data/organelle/mito/summary/table
cd ~/data/organelle/mito/summary/table

# manually set orders in `mito_OG.md`
perl -l -MPath::Tiny -e '
    BEGIN {
        @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_}
            map {s/,\w+//; $_}
            map {s/^###\s*//; $_}
            path(q{~/Scripts/withncbi/doc/mito_OG.md})->lines( { chomp => 1});
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
#90 lines, 3 fields

# phylum family genus abbr taxon_id
cat ~/data/organelle/mito/summary/GENUS.csv |
    grep -v "^#" |
    perl -nla -F"," -e 'print join qq{\t}, ($F[8], $F[5], $F[4], $F[9], $F[0], )' |
    sort |
    uniq \
    > abbr.tmp
cat abbr.tmp | datamash check
#91 lines, 5 fields

tsv-join \
    abbr.tmp \
    --data-fields 4 \
    -f length.tmp \
    --key-fields 1 \
    --append-fields 2,3 \
    > list.tmp
cat list.tmp | datamash check
#90 lines, 7 fields

# sort as orders in mito_OG.md
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
#91 lines, 7 fields

rm *.tmp

```

## Statistics of genome alignments

Some genera will be filtered out here.

Criteria:

* Coverage >= 0.5
* Total number of indels >= 100
* Genome D < 0.05

```bash
cd ~/data/organelle/mito/summary/xlsx

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

cat ~/data/organelle/mito/summary/table/genus_all.lst |
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

cp -f Table_alignment_all.xlsx ~/data/organelle/mito/summary/table
cp -f Table_alignment_all.csv ~/data/organelle/mito/summary/table

```

```bash
cd ~/data/organelle/mito/summary/table

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

rm ~/data/organelle/mito/summary/table/Table_alignment_all.[0-9].csv
rm ~/data/organelle/mito/summary/table/group_*csv

#
cd ~/data/organelle/mito/summary/xlsx
cat ~/data/organelle/mito/summary/table/genus.lst |
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

cp -f ~/data/organelle/mito/summary/xlsx/Table_alignment.xlsx ~/data/organelle/mito/summary/table
cp -f ~/data/organelle/mito/summary/xlsx/Table_alignment.csv ~/data/organelle/mito/summary/table

```

## Groups

NCBI Taxonomy tree

```bash
mkdir -p ~/data/organelle/mito/summary/group
cd ~/data/organelle/mito/summary/group

cat ~/data/organelle/mito/summary/table/genus.lst |
    grep -v "^#" |
    perl -e '
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

bash genera_tree.sh > genera.newick

nw_display -w 600 -s genera.newick |
    rsvg-convert -f pdf -o genera.pdf

```

## Phylogenic trees of each genus with outgroup

```bash
mkdir -p ~/data/organelle/mito/summary/trees
cd ~/data/organelle/plastid/summary/trees

cat ~/Scripts/withncbi/doc/mito_OG.md |
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
cd ~/data/organelle/mito/summary/xlsx

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
mkdir -p ~/data/organelle/mito/summary/fig

cd ~/data/organelle/mito/summary/xlsx

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

cat ../table/genus.lst |
    TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        push @data, { name => $_, };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.03,
                y_max2 => 0.03,
                postfix => q{group_1}, })
                or die Template->error;
        }
    ' \
    > cmd_chart_genus.sh

bash cmd_chart_genus.sh

rm ../xlsx/*.csv
mv ../xlsx/*.pdf ../fig

```

