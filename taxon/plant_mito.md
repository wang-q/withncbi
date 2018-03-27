# Process plant mitochondrion genomes

[TOC levels=1-3]: # " "
- [Process plant mitochondrion genomes](#process-plant-mitochondrion-genomes)
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
- [Add lineage information](#add-lineage-information)
    - [Can't get clear taxon information](#cant-get-clear-taxon-information)
- [Filtering based on valid families and genera](#filtering-based-on-valid-families-and-genera)
- [Find a way to name these](#find-a-way-to-name-these)
- [Download sequences and regenerate lineage information.](#download-sequences-and-regenerate-lineage-information)
- [Create alignment plans](#create-alignment-plans)
- [Aligning](#aligning)
    - [Batch running for groups](#batch-running-for-groups)
    - [Self alignments.](#self-alignments)
    - [Alignments of families for outgroups.](#alignments-of-families-for-outgroups)


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

# Find a way to name these

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


Create `mitochondrion_OG.md` for picking outgroups.

Manually edit it then move to `~/Scripts/withncbi/doc/mitochondrion_OG.md`.

```bash
cd ~/data/organelle/mitochondrion_summary

cat mitochondrion.GENUS.csv |
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
    > mitochondrion_OG.md
```


Create alignments without/with outgroups.

```bash
cd ~/data/organelle/mitochondrion_summary

# tab-separated
# name  t   qs
cat mitochondrion.GENUS.csv |
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

cat mitochondrion.ABBR.csv |
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

# name  t   qs  o
cat genus.tsv |
    perl -nl -a -F"\t" -MPath::Tiny -e '
        BEGIN{
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{~/Scripts/withncbi/doc/mitochondrion_OG.md})->lines({ chomp => 1});
            for (@ls) {
                @fs = split(/,/);
                $h{$fs[0]}= $fs[1];
            }
        }

        if (exists $h{$F[0]}) {
            printf qq{%s\t%s\t%s\t%s\n}, $F[0] . q{_OG}, $F[1], $F[2], $h{$F[0]};
        }' \
    > genus_OG.tsv

# every genera
echo -e "mkdir -p ~/data/organelle/mitochondrion.working \ncd ~/data/organelle/mitochondrion.working\n" > ../mitochondrion.cmd.txt
cat genus.tsv |
    perl ~/Scripts/withncbi/taxon/cmd_template.pl \
        --seq_dir ~/data/organelle/mitochondrion_genomes \
        --csv_taxon ~/data/organelle/mitochondrion_genomes/mitochondrion_ncbi.csv \
        --parallel 8 \
    >> ../mitochondrion.cmd.txt

echo -e "mkdir -p ~/data/organelle/mitochondrion.working \ncd ~/data/organelle/mitochondrion.working\n" > ../mitochondrion.redo.cmd.txt
cat genus.tsv |
    perl ~/Scripts/withncbi/taxon/cmd_template.pl \
        --csv_taxon ~/data/organelle/mitochondrion_genomes/mitochondrion_ncbi.csv \
        --parallel 8 \
    >> ../mitochondrion.redo.cmd.txt

# this is for finding outgroups
echo -e "mkdir -p ~/data/organelle/mitochondrion_families \ncd ~/data/organelle/mitochondrion_families\n" > ../mitochondrion_families.cmd.txt
cat family.tsv |
    perl -n -e '/,\w+,/ and print' |
    perl ~/Scripts/withncbi/taxon/cmd_template.pl \
        --seq_dir ~/data/organelle/mitochondrion_genomes \
        --csv_taxon ~/data/organelle/mitochondrion_genomes/mitochondrion_ncbi.csv \
        --parallel 8 \
    >> ../mitochondrion_families.cmd.txt

# genera with outgroups
echo -e "mkdir -p ~/data/organelle/mitochondrion_OG \ncd ~/data/organelle/mitochondrion_OG\n" > ../mitochondrion_OG.cmd.txt
cat genus_OG.tsv |
    perl ~/Scripts/withncbi/taxon/cmd_template.pl \
        --seq_dir ~/data/organelle/mitochondrion_genomes \
        --csv_taxon ~/data/organelle/mitochondrion_genomes/mitochondrion_ncbi.csv \
        --parallel 8 \
    >> ../mitochondrion_OG.cmd.txt
```

# Aligning

## Batch running for groups

```bash
mkdir -p ~/data/organelle/mitochondrion.working
cd ~/data/organelle/mitochondrion.working

bash ../mitochondrion.cmd.txt 2>&1 | tee log_cmd.txt
# bash ../mitochondrion.redo.cmd.txt 2>&1 | tee log_redo_cmd.txt # skip real_chr and repeatmasker

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

## Self alignments.

```bash
cd ~/data/organelle/

perl -p -e '
    s/mitochondrion\.working/mitochondrion_self.working/g;
    s/multi_batch/self_batch/g;
    s/(\-\-parallel)/--length 1000 \1/g;
' mitochondrion.cmd.txt > mitochondrion_self.cmd.txt

mkdir -p ~/data/organelle/mitochondrion_self.working
cd ~/data/organelle/mitochondrion_self.working

time bash ../mitochondrion_self.cmd.txt 2>&1 | tee log_cmd.txt

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

## Alignments of families for outgroups.

```bash
mkdir -p ~/data/organelle/mitochondrion_families
cd ~/data/organelle/mitochondrion_families

time bash ../mitochondrion_families.cmd.txt 2>&1 | tee log_cmd.txt

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

find ~/data/organelle/mitochondrion_families -type f -path "*_phylo*" -name "*.nwk"

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel -r rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel -r rm

```

