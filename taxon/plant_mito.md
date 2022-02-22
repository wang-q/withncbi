# Plant mitochondrion genomes

- [Plant mitochondrion genomes](#plant-mitochondrion-genomes)
- [Preparation](#preparation)
  * [Download software](#download-software)
  * [Update taxonomy database](#update-taxonomy-database)
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
  * [Restrict taxonomy ids to green plants](#restrict-taxonomy-ids-to-green-plants)
- [Add lineage information](#add-lineage-information)
- [Filtering based on valid families and genera](#filtering-based-on-valid-families-and-genera)
- [Find a way to name these](#find-a-way-to-name-these)
  * [Numbers for higher ranks](#numbers-for-higher-ranks)
- [Download sequences and regenerate lineage information](#download-sequences-and-regenerate-lineage-information)
  * [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
- [Prepare sequences for lastz](#prepare-sequences-for-lastz)
- [Aligning without outgroups](#aligning-without-outgroups)
  * [Create `mito_t_o.md` for picking targets and outgroups.](#create-mito_t_omd-for-picking-targets-and-outgroups)
  * [Create alignments plans without outgroups](#create-alignments-plans-without-outgroups)
  * [Plans for align-able targets](#plans-for-align-able-targets)
  * [Aligning w/o outgroups](#aligning-wo-outgroups)
  * [Aligning with outgroups](#aligning-with-outgroups)
  * [Self alignments](#self-alignments)
- [Summary](#summary)
  * [Copy xlsx files](#copy-xlsx-files)
  * [Genome list](#genome-list)
  * [Statistics of genome alignments](#statistics-of-genome-alignments)
  * [Groups](#groups)
  * [Phylogenic trees of each genus with outgroup](#phylogenic-trees-of-each-genus-with-outgroup)
  * [d1, d2](#d1-d2)

# Preparation

## Download software

```shell
brew install miller librsvg
brew install mash newick_utils
brew install wang-q/tap/nwr wang-q/tap/tsv-utils

```

Install [egaz](https://github.com/wang-q/App-Egaz#installation)

## Update taxonomy database

```shell
rm -fr ~/.nwr

nwr download

nwr txdb

```

*[Update](../db/README.md#get-data-from-ncbi) `~/data/NCBI/taxdmp` before running `strain_info.pl`.*

# Scrap id and acc from NCBI

Use `taxon/gb_taxon_locus.pl` to extract information from RefSeq files.

```shell
mkdir -p ~/data/mito/GENOMES
cd ~/data/mito/GENOMES

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz

gzip -dcf mitochondrion.*.genomic.gbff.gz > genomic.gbff

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv

rm genomic.gbff
#rm *.genomic.gbff.gz

cat refseq_id_seq.csv | grep -v "^#" | wc -l
# 12290

# combine
cat refseq_id_seq.csv |
    sort -u | # duplicated id-seq pair
    sort -t, -k1,1 |
    mlr --icsv --otsv cat \
    > id_seq.tsv

cat id_seq.tsv | grep -v "^#" | wc -l
# 12290

```

## Restrict taxonomy ids to green plants

```text
Eukaryota (2759)
    Viridiplantae (33090)
        Chlorophyta (3041)
        Streptophyta (35493)
```

```shell
cd ~/data/mito/GENOMES

# Viridiplantae 33090
echo -e '#tax_id\taccession' > plant_id_seq.tsv
cat id_seq.tsv |
    nwr restrict 33090 -f stdin -c 1 \
    >> plant_id_seq.tsv

cat plant_id_seq.tsv | grep -v "^#" | wc -l
# 380

cat plant_id_seq.tsv |
    cut -f 1 |
    sort -n |
    tsv-uniq --number --repeated |
    nwr append stdin
#3659    2       Cucumis sativus
#3659    3       Cucumis sativus
#3708    2       Brassica napus
#51329   2       Polytomella parva
#351366  2       Polytomella piriformis

```

Look inside `plant_id_seq.tsv` and remove redundancies

```shell
# Cucumis sativus has 3 chromosomes
# Brassica napus linear plasmid NC_004946
# Polytomella parva has 2 chromosomes
# Polytomella sp. SAG 63-10 has 2 chromosomes

sed -i".bak" "/,NC_004946$/d" plant_id_seq.csv # Brassica napus

```

# Add lineage information

Give ids better shapes for manually checking and automatic filtering.

If you sure, you can add or delete lines and contents in `CHECKME.tsv`.

```shell
mkdir -p ~/data/mito/summary
cd ~/data/mito/summary

# generate a TSV file for manually checking
cat ../GENOMES/plant_id_seq.tsv |
    nwr append stdin |
    nwr append stdin -r species -r genus -r family -r order -r class -r phylum |
    keep-header -- sort -k9,9 -k8,8 -k7,7 -k6,6 -k5,5 \
    > CHECKME.tsv

```

Manually correct lineages.

Split Streptophyta according to classical plant classification.

* Streptophyta
  * Streptophytina
    * Embryophyta
      * Tracheophyta
        * Euphyllophyta
          * Spermatophyta
            * Magnoliopsida 3398 - Angiosperm
            * Acrogymnospermae 1437180 - Gymnosperm
          * Polypodiopsida 241806 - ferns
        * Lycopodiopsida 1521260 - clubmosses
      * Anthocerotophyta 13809 - hornworts
      * Bryophyta 3208 - mosses
      * Marchantiophyta 3195 - liverworts

```shell
cd ~/data/mito/summary

#perl -Mojo -e '
#    g(q{http://www.theplantlist.org/browse/A/})->dom
#    ->find(q{li > a > i[class=family]})
#    ->each( sub { print shift->text . "\n" } );
#    ' > Angiosperms.tmp

# darwin (bsd) need "" for -i
sed -i".bak" CHECKME.csv

# Angiosperms
nwr member 3398 -r family |
    tsv-select -f 2 |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({}\t\w+\t\w+)\tStreptophyta/\1\tAngiosperms/g
        '\'' CHECKME.tsv
    '

# Gymnosperms
nwr member 1437180 -r family |
    tsv-select -f 2 |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({}\t\w+\t\w+)\tStreptophyta/\1\tGymnosperms/g
        '\'' CHECKME.tsv
    '

# Pteridophytes
nwr member 241806 -r family |
    tsv-select -f 2 |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({}\t\w+\t\w+)\tStreptophyta/\1\tPteridophytes/g
        '\'' CHECKME.tsv
    '
nwr member 1521260 -r family |
    tsv-select -f 2 |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({}\t\w+\t\w+)\tStreptophyta/\1\tPteridophytes/g
        '\'' CHECKME.tsv
    '

# Bryophytes
nwr member 13809 -r family |
    tsv-select -f 2 |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({}\t\w+\t\w+)\tStreptophyta/\1\tBryophytes/g
        '\'' CHECKME.tsv
    '

nwr member 3208 -r family |
    tsv-select -f 2 |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({}\t\w+\t\w+)\tStreptophyta/\1\tBryophytes/g
        '\'' CHECKME.tsv
    '

nwr member 3195 -r family |
    tsv-select -f 2 |
    parallel -r -j 1 '
        perl -pi -e '\''
            s/({}\t\w+\t\w+)\tStreptophyta/\1\tBryophytes/g
        '\'' CHECKME.tsv
    '

# Charophyta 轮藻
parallel -j 1 '
    sed -i".bak" "s/{}\tStreptophyta/{}\tCharophyta/" CHECKME.tsv
    ' ::: \
        Charophyceae Chlorokybophyceae Coleochaetophyceae \
        Zygnemophyceae

# Chlorophyta 绿藻
parallel -j 1 '
    sed -i".bak" "s/{}\tStreptophyta/{}\tChlorophyta/" CHECKME.tsv
    ' ::: \
        Mesostigmatophyceae Klebsormidiophyceae

# Ochrophyta 褐藻
parallel -j 1 '
    sed -i".bak" "s/{}\tNA/{}\tOchrophyta/" CHECKME.tsv
    ' ::: \
        Bolidophyceae Dictyochophyceae Eustigmatophyceae \
        Pelagophyceae Phaeophyceae Raphidophyceae \
        Synurophyceae Xanthophyceae

# Rhodophyta 红藻
parallel -j 1 '
    sed -i".bak" "s/{}\tNA/{}\tRhodophyta/" CHECKME.tsv
    ' ::: \
        Bangiophyceae Compsopogonophyceae Florideophyceae \
        Rhodellophyceae Stylonematophyceae

# Cryptophyta 隐藻
parallel -j 1 '
    sed -i".bak" "s/{}\tNA/{}\tCryptophyta/" CHECKME.tsv
    ' ::: \
        Cryptophyta Cryptophyceae

# missing phylums
sed -i".bak" "s/Glaucocystophyceae\tNA/Glaucocystophyceae\tGlaucophyta/" CHECKME.tsv
sed -i".bak" "s/Dinophyceae\tNA/Dinophyceae\tDinoflagellata/" CHECKME.tsv

rm *.tmp *.bak

```

# Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

```text
380 ---------> 372 ---------> 160 ---------> 210
        NA           genus         family
```

```shell
mkdir -p ~/data/mito/summary
cd ~/data/mito/summary

cat CHECKME.tsv | grep -v "^#" | wc -l
# 380

# filter out accessions without linage information (strain, species, genus and family)
cat CHECKME.tsv |
    perl -nla -F"\t" -e '
        /^#/ and next;
        ($F[2] eq q{NA} or $F[3] eq q{NA} or $F[4] eq q{NA} or $F[5] eq q{NA} ) and next;
        print
    ' \
    > valid.tmp

wc -l valid.tmp
# 372

#----------------------------#
# Genus
#----------------------------#
# valid genera
cat valid.tmp |
    perl -nla -F"\t" -e '
        $seen{$F[4]}++;
        END {
            for $k (sort keys %seen) {
                printf qq{\t%s\t\n}, $k if $seen{$k} > 1
            }
        }
    ' \
    > genus.tmp

# intersect between two files
grep -F -f genus.tmp valid.tmp > valid.genus.tmp

wc -l valid.genus.tmp
# 160

#----------------------------#
# Family
#----------------------------#
# get some genera back as candidates for outgroup
cat valid.genus.tmp |
    perl -nla -F"\t" -e 'printf qq{\t$F[5]\t\n}' \
    > family.tmp

# intersect between two files
grep -F -f family.tmp valid.tmp > valid.family.tmp

wc -l valid.family.tmp
# 210

#----------------------------#
# results produced in this step
#----------------------------#
head -n 1 CHECKME.tsv > DOWNLOAD.tsv
cat valid.family.tmp |
    sort -t$'\t' -k9,9 -k8,8 -k7,7 -k6,6 -k5,5 \
    >> DOWNLOAD.tsv

# clean
rm *.tmp *.bak

```

# Find a way to name these

Seems it's OK to use species as names.

```shell
cd ~/data/mito/summary

# sub-species
cat DOWNLOAD.tsv |
    perl -nl -a -F"\t" -e '
        /^#/i and next;
        $seen{$F[3]}++;
        END {
            for $k (keys %seen){printf qq{%s\t%d\n}, $k, $seen{$k} if $seen{$k} > 1}
        };
    ' |
    sort
#Beta vulgaris   2
#Brassica napus  2
#Brassica rapa   2
#Cucumis sativus 3
#Oryza sativa    2
#Polytomella parva       2
#Polytomella piriformis  2
#Zea mays        2

# strain name not equal to species
cat DOWNLOAD.tsv |
    grep -v '^#' |
    tsv-filter --ff-str-ne 3:4 |
    tsv-select -f 3 |
    sort
#Aegilops speltoides var. ligustica
#Beta vulgaris subsp. maritima
#Beta vulgaris subsp. vulgaris
#Brassica rapa subsp. oleifera
#Calypogeia fissa subsp. neogaea
#Marchantia polymorpha subsp. ruderalis
#Oryza sativa Indica Group
#Oryza sativa Japonica Group
#Zea mays subsp. mays
#Zea mays subsp. parviglumis

```

Create abbreviations.

```shell
cd ~/data/mito/summary

head -n 1 DOWNLOAD.tsv |
    sed 's/$/\tabbr/' \
     > ABBR.tsv

cat DOWNLOAD.tsv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s '\t' -m 0 --shortsub |
    mlr --icsv --otsv cat |
    sort -t$'\t' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> ABBR.tsv

```

## Numbers for higher ranks

21 orders, 22 families, 36 genera and 106 species.

```shell
cd ~/data/mito/summary

## valid genera
#cat ABBR.tsv |
#    grep -v "^#" |
#    perl -nl -a -F"\t" -e '
#        $seen{$F[4]}++;
#        END {
#            for $k (sort keys %seen) {
#                printf qq{\t%s\t\n}, $k if $seen{$k} > 1
#            }
#        }
#    ' \
#    > genus.valid.tmp

cat ABBR.tsv |
    tsv-uniq -H -f genus --at-least 2 |
    tsv-select -H -f genus |
    sed 's/^/\t/' |
    sed 's/$/\t/' \
    > genus.valid.tmp

# intersect between two files
grep -F -f genus.valid.tmp ABBR.tsv > GENUS.tsv

wc -l GENUS.tsv
# 161

# count every ranks
for rank in species genus family order; do
    cat GENUS.tsv |
        tsv-summarize -H --group-by ${rank} --count |
        keep-header -- wc -l
done
#species count
#151
#genus   count
#45
#family  count
#28
#order   count
#24

# sort by multiply columns: phylum, order, family, genus, accession
# Older accessions have better sequencing qualities
cat GENUS.tsv |
    keep-header -- sort -t$'\t' -k9,9 -k7,7 -k6,6 -k5,5 -k2,2 \
    >> GENUS.tmp
mv GENUS.tmp GENUS.tsv

# clean
rm *.tmp *.bak

```

# Download sequences and regenerate lineage information

```shell
cd ~/data/mito/GENOMES

echo "#strain_name,accession,strain_taxon_id" > name_acc_id.csv
cat ../summary/ABBR.tsv |
    grep -v '^#' |
    perl -nl -a -F"\t" -e 'print qq{$F[9],$F[1],$F[0]}' |
    sort \
    >> name_acc_id.csv

# Local, Runtime 4 seconds.
cat ../summary/ABBR.tsv |
    grep -v '^#' |
    perl -nla -F"\t" -e 'print qq{$F[0],$F[9]}' |
    uniq |
    perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file taxon_ncbi.csv

# Some warnings about trans-splicing genes from BioPerl, just ignore them
# eutils restricts 3 requests per second
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

## Raw phylogenetic tree by MinHash

```shell
mkdir -p ~/data/mito/mash
cd ~/data/mito/mash

for name in $(tsv-select ../summary/ABBR.tsv -H -f abbr | sed -e '1d'); do
    2>&1 echo "==> ${name}"

    if [[ -e ${name}.msh ]]; then
        continue
    fi

    find ../GENOMES/${name} -name "*.fa" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 4 - -I "${name}" -o ${name}
done

```

* h=0.1 as we need outgroup ~ 0.05

```shell
cd ~/data/mito/summary
mash triangle -E -p 8 -l <(
        tsv-select ABBR.tsv -H -f abbr |
            sed '1d' |
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

        group <- cutree(clusters, h=0.1)
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

nw_display -s -b 'visibility:hidden' -w 600 -v 30 tree.nwk |
    rsvg-convert -o plant_mito.png

```

* Abnormal strains

```shell
cd ~/data/mito/summary

# genus
tsv-select GENUS.tsv -H -f genus | sed -e '1d' | uniq |
    parallel -j 1 -k '
        group=$(
            tsv-join groups.tsv -d 2 \
                -f <(cat GENUS.tsv | grep -w {} | tsv-select -f 10) \
                -k 1 |
                cut -f 1 |
                sort |
                uniq
        )
        number=$(echo "${group}" | wc -l)
        echo -e "{}\t${number}"
    ' |
    tsv-join -d 1 -f GENUS.tsv -k 5 -a 9 |
    tsv-filter --ne 2:1
#Solanum 2       Angiosperms
#Bracteacoccus   2       Chlorophyta
#Caulerpa        2       Chlorophyta
#Chlamydomonas   3       Chlorophyta
#Polytomella     4       Chlorophyta
#Chlorella       3       Chlorophyta
#Prototheca      2       Chlorophyta
#Chloroparvula   2       Chlorophyta
#Chloropicon     5       Chlorophyta
#Dunaliella      2       Chlorophyta
#Ulva    7       Chlorophyta

# family
tsv-select GENUS.tsv -H -f family | sed -e '1d' | uniq |
    parallel -j 1 -k '
        group=$(
            tsv-join groups.tsv -d 2 \
                -f <(cat GENUS.tsv | grep -w {} | tsv-select -f 10) \
                -k 1 |
                cut -f 1 |
                sort |
                uniq
        )
        number=$(echo "${group}" | wc -l)
        echo -e "{}\t${number}"
    ' |
    tsv-join -d 1 -f GENUS.tsv -k 6 -a 9 |
    tsv-filter --ne 2:1
#Fabaceae        4       Angiosperms
#Malvaceae       2       Angiosperms
#Poaceae 3       Angiosperms
#Solanaceae      2       Angiosperms
#Bracteacoccaceae        2       Chlorophyta
#Caulerpaceae    2       Chlorophyta
#Chlamydomonadaceae      7       Chlorophyta
#Chlorellaceae   5       Chlorophyta
#Chloropicaceae  7       Chlorophyta
#Dunaliellaceae  2       Chlorophyta
#Ulvaceae        7       Chlorophyta

```

# Prepare sequences for lastz

```shell
cd ~/data/mito/GENOMES

find . -maxdepth 1 -mindepth 1 -type d |
    sort |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo >&2 "==> {}"

        if [ -e {}/chr.fasta ]; then
            echo >&2 "    {} has been processed";
            exit;
        fi

        egaz prepseq \
            {} \
            --gi -v --repeatmasker " --gff --parallel 4"
    '

# restore to original states
#for suffix in .2bit .fasta .fasta.fai .sizes .rm.out .rm.gff; do
#    find . -name "*${suffix}" | parallel --no-run-if-empty rm
#done

```

# Aligning without outgroups

## Create `mito_t_o.md` for picking targets and outgroups.

Manually edit it then move to `~/Scripts/withncbi/doc/mito_t_o.md`.

* Listed targets were well curated.

* Outgroups can be changes with less intentions.

```shell
cd ~/data/mito/summary

cat GENUS.tsv |
    grep -v "^#" |
    perl -na -F"\t" -e '
        BEGIN{
            our ($phylum, $family, $genus) = (q{}, q{}, q{});
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

# 45
cat mito_t_o.md |
    grep -v '^#' |
    grep -E '\S+' |
    wc -l

# mv mito_t_o.md ~/Scripts/withncbi/doc/mito_t_o.md

```

## Create alignments plans without outgroups

```text
GENUS.tsv
#strain_taxon_id accession strain species genus family order class phylum abbr
```

**34** genera and **19** families.

```shell
mkdir -p ~/data/mito/taxon
cd ~/data/mito/taxon

echo -e "#Serial\tGroup\tCount\tTarget" > group_target.tsv

# genus
cat ../summary/GENUS.tsv |
    grep -v "^#" |
    SERIAL=1 perl -na -F"\t" -MPath::Tiny -e '
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

# family
cat ../summary/ABBR.tsv |
    grep -v "^#" |
    SERIAL=1001 perl -na -F"\t" -MPath::Tiny -e '
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

```shell
cd ~/data/mito/taxon

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
    SERIAL=2001 perl -na -F"\t" -MPath::Tiny -e '
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

Start a MySQL server, which aligndb needs. `mysqld_safe`

```shell
cd ~/data/mito/

# genus
cat taxon/group_target.tsv |
    tsv-filter -H  --ge 1:1 --le 1:1000 |
    sed -e '1d' | #grep -w "^24" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/genus/{2} \
            --taxon ~/data/mito/GENOMES/taxon_ncbi.csv \
            --rawphylo --aligndb --parallel 4 -v

        bash groups/genus/{2}/1_pair.sh
        bash groups/genus/{2}/2_rawphylo.sh
        bash groups/genus/{2}/3_multi.sh
        bash groups/genus/{2}/6_chr_length.sh
        bash groups/genus/{2}/7_multi_aligndb.sh
    '

# family
cat taxon/group_target.tsv |
    tsv-filter -H --ge 1:1001 --le 1:2000 |
    sed -e '1d' | #grep -w "^1008" |
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

# mash
cat taxon/group_target.tsv |
    tsv-filter -H --ge 1:2001 --le 1:3000 |
    sed -e '1d' | #grep -w "^2001" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/mash/{2} \
            --rawphylo --parallel 4 -v

        bash groups/mash/{2}/1_pair.sh
        bash groups/mash/{2}/2_rawphylo.sh
        bash groups/mash/{2}/3_multi.sh
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

* Abnormal groups

```shell
cd ~/data/mito/

find groups/ -name "pairwise.coverage.csv" | sort |
    parallel -j 4 -k '
        cover=$(cat {} | grep -w "intersect" | cut -d, -f 4)

        if [[ "$cover" == "" ]]; then
            cover=$(cat {} | cut -d, -f 4 | tsv-summarize -H --mean 1 | sed -e "1d")
        fi

        taxon=$(
            echo {//} |
                sed -e "s/^groups\///g" -e "s/\/Results//g" |
                sed -e "s/^family\///g" -e "s/^genus\///g" -e "s/^mash\///g"
            )
        phylum=$(
            cat summary/ABBR.tsv |
                grep -w ${taxon} |
                tsv-select -f 9 |
                head -n 1
            )
        group=$(
            echo {//} |
                sed -e "s/^groups\///g" -e "s/\/Results//g" |
                sed -E "s/\/.+//g"
            )

        if [ $(bc <<< "${cover} < 0.5") -eq 1 ]; then
            echo -e "${taxon}\t${cover}\t${phylum}\t${group}"
        fi
    '

#Asteraceae      0.1488  Angiosperms     family
#Brassicaceae    0.1181  Angiosperms     family
#Chenopodiaceae  0.1791  Angiosperms     family
#Chlorellaceae   0.09496 Chlorophyta     family
#Chloropicaceae  0.0317  Chlorophyta     family
#Cucurbitaceae   0.0140  Angiosperms     family
#Fabaceae        0.0000  Angiosperms     family
#Malvaceae       0.0103  Angiosperms     family
#Nyctaginaceae   0.3327  Angiosperms     family
#Poaceae 0.0000  Angiosperms     family
#Rutaceae        0.0124  Angiosperms     family
#Salicaceae      0.4175  Angiosperms     family
#Solanaceae      0.0000  Angiosperms     family
#Ulvaceae        0.0744  Chlorophyta     family
#Bracteacoccus   0.2154  Chlorophyta     genus
#Caulerpa        0.0139  Chlorophyta     genus
#Chlorella       0.1734  Chlorophyta     genus
#Chloroparvula   0.1954  Chlorophyta     genus
#Chloropicon     0.0940  Chlorophyta     genus
#Corchorus       0.0548  Angiosperms     genus
#Prototheca      0.0866  Chlorophyta     genus
#Senna   0.3345  Angiosperms     genus
#Solanum 0.0371  Angiosperms     genus
#Ulva    0.0811  Chlorophyta     genus
#group_23        0.3832          mash

```

## Aligning with outgroups

* Review alignments and phylogenetic trees generated in `groups/family/` and `groups/mash/`

* Add outgroups to `mito_t_o.md` manually.

* *D* between target and outgroup should be around **0.05**.

```shell
cd ~/data/mito/

# genus_og
cat taxon/group_target.tsv |
    tsv-filter -H --le 1:1000 |
    sed -e '1d' | #grep -w "^24" |
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
            --multi -o groups/genus_og/{2}_og \
            --outgroup ${outgroup} \
            --taxon ~/data/mito/GENOMES/taxon_ncbi.csv \
            --rawphylo --aligndb --parallel 4 -v

        bash groups/genus_og/{2}_og/1_pair.sh
        bash groups/genus_og/{2}_og/2_rawphylo.sh
        bash groups/genus_og/{2}_og/3_multi.sh
        bash groups/genus_og/{2}_og/6_chr_length.sh
        bash groups/genus_og/{2}_og/7_multi_aligndb.sh
    '

```

## Self alignments

```shell
cd ~/data/mito/

cat taxon/group_target.tsv |
    tsv-filter -H --le 1:1000 |
    sed -e '1d' | # grep -w "^24" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --self -o groups/self/{2} \
            --circos --parallel 4 -v

        bash groups/self/{2}/1_self.sh
        bash groups/self/{2}/3_proc.sh
        bash groups/self/{2}/4_circos.sh
    '

```

# Summary

## Copy xlsx files

```shell
mkdir -p ~/data/mito/summary/xlsx
cd ~/data/mito/summary/xlsx

find ../../groups/genus -type f -name "*.common.xlsx" |
    grep -v "vs[A-Z]" |
    parallel 'cp {} .'

find ../../groups/genus_og -type f -name "*.common.xlsx" |
    grep -v "vs[A-Z]" |
    parallel 'cp {} .'

```

## Genome list

Create `list.csv` from `GENUS.tsv` with sequence lengths.

```shell
mkdir -p ~/data/mito/summary/table
cd ~/data/mito/summary/table

# manually set orders in `mito_t_o.md`
perl -l -MPath::Tiny -e '
    BEGIN {
        @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_}
            map {s/,.+$//; $_}
            map {s/^###\s*//; $_}
            path(q{~/Scripts/withncbi/doc/mito_t_o.md})->lines({chomp => 1});
    }
    for (@ls) {
        (/^\s*$/ or /^##\s+/ or /^#\s+(\w+)/) and next;
        print $_
    }
    ' \
    > genus_all.lst

# abbr accession length
find ../../groups/genus -type f -name "chr_length.csv" |
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
#123 lines, 3 fields

# phylum family genus abbr taxon_id
cat ~/data/mito/summary/GENUS.tsv |
    grep -v "^#" |
    perl -nla -F"\t" -e 'print join qq{\t}, ($F[8], $F[5], $F[4], $F[9], $F[0], )' |
    sort |
    uniq \
    > abbr.tmp
cat abbr.tmp | datamash check
#155 lines, 5 fields

tsv-join \
    abbr.tmp \
    --data-fields 4 \
    -f length.tmp \
    --key-fields 1 \
    --append-fields 2,3 \
    > list.tmp
cat list.tmp | datamash check
#123 lines, 7 fields

# sort as orders in mito_t_o.md
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
#124 lines, 7 fields

rm *.tmp

```

## Statistics of genome alignments

Some genera will be filtered out here.

Criteria:

* Coverage >= 0.5
* Total number of indels >= 100
* D of multiple alignments < 0.05

```shell
cd ~/data/mito/summary/xlsx

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

cat ~/data/mito/summary/table/genus_all.lst |
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

cp -f Table_alignment_all.xlsx ~/data/mito/summary/table
cp -f Table_alignment_all.csv ~/data/mito/summary/table

```

```shell
cd ~/data/mito/summary/table

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

rm ~/data/mito/summary/table/Table_alignment_all.[0-9].csv
rm ~/data/mito/summary/table/group_*csv

#
cd ~/data/mito/summary/xlsx
cat ~/data/mito/summary/table/genus.lst |
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

cp -f ~/data/mito/summary/xlsx/Table_alignment.xlsx ~/data/mito/summary/table
cp -f ~/data/mito/summary/xlsx/Table_alignment.csv ~/data/mito/summary/table

```

## Groups

NCBI Taxonomy tree

```shell
mkdir -p ~/data/mito/summary/group
cd ~/data/mito/summary/group

cat ~/data/mito/summary/table/genus.lst |
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
    > genera_tree.sh

bash genera_tree.sh > genera.newick

nw_display -s -b 'visibility:hidden' -w 600 -v 30 genera.newick |
    rsvg-convert -o genera.png

```

## Phylogenic trees of each genus with outgroup

```shell
mkdir -p ~/data/mito/summary/trees
cd ~/data/mito/summary/trees

cat ~/Scripts/withncbi/doc/mito_t_o.md |
    grep -v "^#" |
    grep . |
    cut -d',' -f 1 \
    > list.txt

find ../../groups/genus_og -type f -path "*Results*" -name "*.nwk" |
    grep -v ".raw." |
    parallel -j 1 'cp {} .'

```

## d1, d2

```shell
mkdir -p ~/data/mito/summary/fig
cd ~/data/mito/summary/xlsx

# d1
cat ../table/genus.lst |
    grep -v "^#" |
    parallel -j 1 -k '
        >&2 echo "==> {}"
        plotr xlsx {}.common.xlsx -s d1_pi_gc_cv -o stdout |
            sed "1d" |
            tsv-select -f 1,2 |
            tsv-filter --le 1:5 --ge 1:0 |
            perl -nla -e '\''
                $F[1] = sprintf q{%.6f}, $F[1];
                print join qq{\t}, @F, {};
            '\''
    ' |
    (echo -e "X\tY\tgroup" && cat) \
    > d1.tsv

cat d1.tsv |
    tsv-summarize -H --group-by 1 --mean 2:Y |
    perl -pe 's/$/\tgroup/' \
    > d1.mean.tsv

plotr lines d1.tsv \
    --ymm 0,0.02 \
    --xl "Distance to indels ({italic(d)[1]})" \
    --yl "Nucleotide divergence ({italic(D)})"

plotr lines d1.mean.tsv --style blue \
    --ymm 0,0.02 \
    --xl "Distance to indels ({italic(d)[1]})" \
    --yl "Nucleotide divergence ({italic(D)})"

# d2
cat ../table/genus.lst |
    grep -v "^#" |
    parallel -j 1 -k '
        >&2 echo "==> {}"
        plotr xlsx {}.common.xlsx -s d2_pi_gc_cv -o stdout |
            sed "1d" |
            tsv-select -f 1,2 |
            tsv-filter --le 1:20 --ge 1:0 |
            perl -nla -e '\''
                $F[1] = sprintf q{%.6f}, $F[1];
                print join qq{\t}, @F, {};
            '\''
    ' |
    (echo -e "X\tY\tgroup" && cat) \
    > d2.tsv

cat d2.tsv |
    tsv-summarize -H --group-by 1 --mean 2:Y |
    perl -pe 's/$/\tgroup/' \
    > d2.mean.tsv

plotr lines d2.tsv \
    --ymm 0,0.03 \
    --xl "Reciprocal of indel density ({italic(d)[2]})" \
    --yl "Nucleotide divergence ({italic(D)})"

plotr lines d2.mean.tsv --style blue \
    --ymm 0,0.03 \
    --xl "Reciprocal of indel density ({italic(d)[2]})" \
    --yl "Nucleotide divergence ({italic(D)})"

mv *.pdf ../fig

```

