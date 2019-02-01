# Process plant mitochondrion genomes

[TOC levels=1-3]: # " "
- [Process plant mitochondrion genomes](#process-plant-mitochondrion-genomes)
- [Update taxdmp](#update-taxdmp)
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
    - [Restrict taxonomy ids to green plants with `taxon/id_restrict.pl`.](#restrict-taxonomy-ids-to-green-plants-with-taxonid_restrictpl)
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
    - [Create `mito_OG.md` for picking outgroups](#create-mito_ogmd-for-picking-outgroups)
    - [Create alignments plans with outgroups](#create-alignments-plans-with-outgroups)
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

Got **9355** accessions.

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

# 9355
cat refseq_id_seq.csv | grep -v "^#" | wc -l

# combine
cat webpage_id_seq.csv refseq_id_seq.csv |
    sort -u | # duplicated id-seq pair
    sort -t, -k1,1 \
    > id_seq.csv

# 9355
cat id_seq.csv | grep -v "^#" | wc -l

```

## Restrict taxonomy ids to green plants with `taxon/id_restrict.pl`.

```bash
cd ~/data/organelle/mito/GENOMES

echo '#strain_taxon_id,accession' > plant_id_seq.csv
cat id_seq.csv |
    grep -v "^#" |
    perl ~/Scripts/withncbi/taxon/id_restrict.pl -s "," -a 33090 \
    >> plant_id_seq.csv

# 242
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

Split Streptophyta according to http://www.theplantlist.org/

```bash
cd ~/data/organelle/mito/summary

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

FIXME

# Filtering based on valid families and genera

Species and genus should not be "NA" and genus has 2 or more members.

```text
238 ---------> 231 ---------> 95 ---------> 126
        NA           genus         family
```


```bash
mkdir -p ~/data/organelle/mito/summary
cd ~/data/organelle/mito/summary

# 237
cat CHECKME.csv | grep -v "^#" | wc -l

# filter out accessions without linage information (strain, species, genus and family)
cat CHECKME.csv |
    perl -nla -F"," -e '
        /^#/ and next;
        ($F[2] eq q{NA} or $F[3] eq q{NA} or $F[4] eq q{NA} or $F[5] eq q{NA} ) and next;
        print
    ' \
    > valid.tmp

# 231
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

# 95
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

# 125
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
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 |
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

Numbers for higher ranks are: 17 orders, 18 families, 31 genera and 88 species.

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

# 95
wc -l GENUS.csv

# count every ranks
#  17 order.list.tmp
#  18 family.list.tmp
#  31 genus.list.tmp
#  88 species.list.tmp
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
cd ~/data/organelle/mito/GENOMES

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
cd ~/data/organelle/mito/summary

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
cd ~/data/organelle/mito/summary

cat <<'EOF' > egaz_template_multi.tt

# [% name %]
egaz template \
    ~/data/organelle/mito/GENOMES/[% t %] \
[% FOREACH q IN qs -%]
    ~/data/organelle/mito/GENOMES/[% q %] \
[% END -%]
[% IF o -%]
    ~/data/organelle/mito/GENOMES/[% o %] \
    --outgroup [% o %] \
[% END -%]
    --multi -o [% name %] \
    --taxon ~/data/organelle/mito/GENOMES/taxon_ncbi.csv \
    --rawphylo --aligndb --parallel 8 -v

EOF

# every genera
echo "mkdir -p ~/data/organelle/mito/genus"  > ../cmd.txt
echo "cd       ~/data/organelle/mito/genus" >> ../cmd.txt
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
echo "mkdir -p ~/data/organelle/mito/family"  > ../family.cmd.txt
echo "cd       ~/data/organelle/mito/family" >> ../family.cmd.txt
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
mkdir -p ~/data/organelle/mito/genus
cd ~/data/organelle/mito/genus

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
mkdir -p ~/data/organelle/mito/family
cd ~/data/organelle/mito/family

time bash ../family.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing ${d} <====\""
    echo bash ${d}/1_pair.sh;
    echo bash ${d}/2_rawphylo.sh;
    echo bash ${d}/3_multi.sh;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

find ~/data/organelle/mito/family -type f -name "*.nwk"

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw"   | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

# Aligning with outgroups

## Create `mito_OG.md` for picking outgroups

Manually edit it then move to `~/Scripts/withncbi/doc/mito_OG.md`.

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
            printf qq{%s\n}, $genus;
        }
    ' \
    > mito_OG.md

```

## Create alignments plans with outgroups

```bash
cd ~/data/organelle/mito/summary

# name  t   qs  o
cat genus.tsv |
    perl -nla -F"\t" -MPath::Tiny -e '
        BEGIN{
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{~/Scripts/withncbi/doc/mito_OG.md})->lines({ chomp => 1});
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
echo "mkdir -p ~/data/organelle/mito/OG"  > ../OG.cmd.txt
echo "cd       ~/data/organelle/mito/OG" >> ../OG.cmd.txt
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

In previous steps, we have manually edited `~/Scripts/withncbi/doc/mito_OG.md` and generated
`genus_OG.tsv`.

*D* between target and outgroup should be around **0.05**.

```bash
mkdir -p ~/data/organelle/mito/OG
cd ~/data/organelle/mito/OG

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

