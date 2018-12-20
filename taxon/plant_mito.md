# Process plant mitochondrion genomes

[TOC levels=1-3]: # " "
- [Process plant mitochondrion genomes](#process-plant-mitochondrion-genomes)
- [Scrap id and acc from NCBI](#scrap-id-and-acc-from-ncbi)
- [Add lineage information](#add-lineage-information)
    - [Can't get clear taxon information](#cant-get-clear-taxon-information)
- [Filtering based on valid families and genera](#filtering-based-on-valid-families-and-genera)
- [Find a way to name these](#find-a-way-to-name-these)
- [Download sequences and regenerate lineage information.](#download-sequences-and-regenerate-lineage-information)
- [Prepare sequences for lastz](#prepare-sequences-for-lastz)
- [Create alignment plans](#create-alignment-plans)
- [Aligning](#aligning)
    - [Batch running for groups](#batch-running-for-groups)
    - [Alignments of families for outgroups.](#alignments-of-families-for-outgroups)
    - [Self alignments](#self-alignments)


# Scrap id and acc from NCBI

Open browser and visit
[NCBI mitochondrion page](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=33090&opt=organelle).
Save page to a local file, html only. In this case, it's
`doc/green_plants_mitochondrion_180325.html`.

All [Eukaryota](https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=2759&opt=organelle),
`doc/eukaryota_mitochondrion_180325.html`.

```text
Eukaryota (2759)                8746
    Viridiplantae (33090)       221
        Chlorophyta (3041)      51
        Streptophyta (35493)    170
```

Use `taxon/id_seq_dom_select.pl` to extract Taxonomy ids and genbank accessions from all history
pages.

Got **217** accessions.

```bash
mkdir -p ~/data/organelle/mito/GENOMES
cd ~/data/organelle/mito/GENOMES

rm webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/eukaryota_mitochondrion_181207.html \
    >> webpage_id_seq.csv
perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/eukaryota_mitochondrion_180325.html \
    >> webpage_id_seq.csv

perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/green_plants_mitochondrion_181207.html \
    >> webpage_id_seq.csv    
perl ~/Scripts/withncbi/taxon/id_seq_dom_select.pl \
    ~/Scripts/withncbi/doc/green_plants_mitochondrion_180325.html \
    >> webpage_id_seq.csv    

```

Use `taxon/gb_taxon_locus.pl` to extract information from refseq genbank files.

```bash
cd ~/data/organelle/mito/GENOMES

wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.1.genomic.gbff.gz
wget -N ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/mitochondrion.2.genomic.gbff.gz

gzip -dcf mitochondrion.*.genomic.gbff.gz > mitochondrion.genomic.gbff

perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl mitochondrion.genomic.gbff > refseq_id_seq.csv

rm mitochondrion.genomic.gbff

# 9103
cat refseq_id_seq.csv | grep -v "^#" | wc -l

# combine
cat webpage_id_seq.csv refseq_id_seq.csv |
    sort -u | # duplicated id-seq pair
    sort -t, -k1,1 \
    > mitochondrion_id_seq.csv

# 9132
cat mitochondrion_id_seq.csv | grep -v "^#" | wc -l

```

Restrict taxonomy ids to green plants with `taxon/id_restrict.pl`.

```bash
cd ~/data/organelle/mito/GENOMES

echo '#strain_taxon_id,accession' > plant_mitochondrion_id_seq.csv
cat mitochondrion_id_seq.csv |
    grep -v "^#" |
    perl ~/Scripts/withncbi/taxon/id_restrict.pl -s "," -a 33090 \
    >> plant_mitochondrion_id_seq.csv

# 232
cat plant_mitochondrion_id_seq.csv | grep -v "^#" | wc -l

cat plant_mitochondrion_id_seq.csv |
    cut -d',' -f 1 |
    sort -n |
    uniq -c |
    grep -v -E '\s+1\s+'
#      3 3659 Cucumis sativus has 3 chromosomes
#      2 3702 Arabidopsis thaliana NC_001284.2 was removed by RefSeq staff
#      2 3708 Brassica napus linear plasmid NC_004946
#      2 39946 Oryza sativa indica plasmid B2 NC_001776
#      2 39947 Oryza sativa japonica plasmid B1 NC_001751
#      2 51329 Polytomella parva has 2 chromosomes
#      2 351366 Polytomella sp. SAG 63-10 has 2 chromosomes

sed -i".bak" "/,NC_001284$/d" plant_mitochondrion_id_seq.csv # Arabidopsis thaliana
sed -i".bak" "/,NC_004946$/d" plant_mitochondrion_id_seq.csv # Brassica napus
sed -i".bak" "/,NC_001751$/d" plant_mitochondrion_id_seq.csv # Oryza sativa japonica
sed -i".bak" "/,NC_001776$/d" plant_mitochondrion_id_seq.csv # Oryza sativa indica

# Vicia faba mitochondrial plasmid MtVFPL3 NC_011084
# Zea mays mitochondrial plasmid pBMSmt1.9 NC_001400
sed -i".bak" "/,NC_011084$/d" plant_mitochondrion_id_seq.csv
sed -i".bak" "/,NC_001400$/d" plant_mitochondrion_id_seq.csv

```

# Add lineage information

Give ids better shapes for manually checking and automatic filtering.

*Update `~/data/NCBI/taxdmp` before running `id_project_to.pl`*.

If you sure, you can add or delete lines and contents in `mitochondrion.CHECKME.csv`.

```bash
mkdir -p ~/data/organelle/mito/summary
cd ~/data/organelle/mito/summary

# generate a .csv file for manually checking
echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum' > mitochondrion.CHECKME.csv
cat ../GENOMES/plant_mitochondrion_id_seq.csv |
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
cd ~/data/organelle/mito/summary

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
224 ---------> 218 ---------> 88 ---------> 116
        NA           genus         family
```


```bash
mkdir -p ~/data/organelle/mito/summary
cd ~/data/organelle/mito/summary

# filter out accessions without linage information (strain, species, genus and family)
cat mitochondrion.CHECKME.csv |
    perl -nla -F"," -e '
        /^#/ and next;
        ($F[2] eq q{NA} or $F[3] eq q{NA} or $F[4] eq q{NA} or $F[5] eq q{NA} ) and next;
        print
    ' \
    > mitochondrion.tmp

# 218
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

# 88
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

# 116
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
cd ~/data/organelle/mito/summary

echo '#strain_taxon_id,accession,strain,species,genus,family,order,class,phylum,abbr' > mitochondrion.ABBR.csv
cat mitochondrion.DOWNLOAD.csv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "3,4,5" -s "," -m 0 |
    sort -t',' -k9,9 -k7,7 -k6,6 -k10,10 \
    >> mitochondrion.ABBR.csv

```

# Download sequences and regenerate lineage information.

```bash
cd ~/data/organelle/mito/GENOMES

echo "#strain_name,accession,strain_taxon_id" > mitochondrion_name_acc_id.csv
cat ../summary/mitochondrion.ABBR.csv |
    grep -v '^#' |
    perl -nl -a -F"," -e 'print qq{$F[9],$F[1],$F[0]}' |
    sort \
    >> mitochondrion_name_acc_id.csv

# local, Runtime 10 seconds.
# with --entrez, Runtime 7 minutes and 23 seconds.
# And which-can't-find is still which-can't-find.
cat ../summary/mitochondrion.ABBR.csv |
    grep -v '^#' |
    perl -nla -F"," -e 'print qq{$F[0],$F[9]}' |
    uniq |
    perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file mitochondrion_ncbi.csv

# some warnings from bioperl, just ignore them
perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -f mitochondrion_name_acc_id.csv \
    2>&1 |
    tee mitochondrion_seq.log

# count downloaded sequences
find . -name "*.fa" | wc -l

```


# Prepare sequences for lastz

```bash
cd ~/data/organelle/mito/GENOMES

find . -maxdepth 1 -type d -path "*/*" |
    sort |
    parallel --no-run-if-empty --linebuffer -k -j 2 '
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

# Create alignment plans

We got **116** accessions.

Numbers for higher ranks are: 16 orders, 17 families, 29 genera and 84 species.

```bash
cd ~/data/organelle/mito/summary

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

# 88
wc -l mitochondrion.GENUS.csv

# count every ranks
#  16 order.list.tmp
#  17 family.list.tmp
#  29 genus.list.tmp
#  84 species.list.tmp
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
cd ~/data/organelle/mito/summary

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
cd ~/data/organelle/mito/summary

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
    perl -nla -F"\t" -MPath::Tiny -e '
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
        }
    ' \
    > genus_OG.tsv

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
    --taxon ~/data/organelle/mito/GENOMES/mitochondrion_ncbi.csv \
    --rawphylo --parallel 8 -v

EOF

# every genera
echo "mkdir -p ~/data/organelle/mito/genus"  > ../mitochondrion.cmd.txt
echo "cd       ~/data/organelle/mito/genus" >> ../mitochondrion.cmd.txt
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
    >> ../mitochondrion.cmd.txt

# this is for finding outgroups
echo "mkdir -p ~/data/organelle/mito/family"  > ../mitochondrion.family.cmd.txt
echo "cd       ~/data/organelle/mito/family" >> ../mitochondrion.family.cmd.txt
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
    >> ../mitochondrion.family.cmd.txt

# genera with outgroups
echo "mkdir -p ~/data/organelle/mito/OG"  > ../mitochondrion.OG.cmd.txt
echo "cd       ~/data/organelle/mito/OG" >> ../mitochondrion.OG.cmd.txt
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
    >> ../mitochondrion.OG.cmd.txt

cat <<'EOF' > egaz_templates_self.tt

# [% name %]
egaz template \
    ~/data/organelle/mito/GENOMES/[% t %] \
[% FOREACH q IN qs -%]
    ~/data/organelle/mito/GENOMES/[% q %] \
[% END -%]
    --self -o [% name %] \
    --taxon ~/data/organelle/mito/GENOMES/mitochondrion_ncbi.csv \
    --circos --aligndb --parallel 8 -v

EOF

# every genera
echo "mkdir -p ~/data/organelle/mito/self"  > ../mitochondrion.self.cmd.txt
echo "cd       ~/data/organelle/mito/self" >> ../mitochondrion.self.cmd.txt
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
    >> ../mitochondrion.self.cmd.txt

```


# Aligning

## Batch running for groups

```bash
mkdir -p ~/data/organelle/mito/genus
cd ~/data/organelle/mito/genus

bash ../mitochondrion.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing ${d} <====\""
    echo bash ${d}/1_pair.sh;
    echo bash ${d}/2_rawphylo.sh;
    echo bash ${d}/3_multi.sh;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

## Alignments of families for outgroups.

```bash
mkdir -p ~/data/organelle/mito/family
cd ~/data/organelle/mito/family

time bash ../mitochondrion.family.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing ${d} <====\""
    echo bash ${d}/1_pair.sh;
    echo bash ${d}/2_rawphylo.sh;
    echo bash ${d}/3_multi.sh;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

find ~/data/organelle/mito/family -type f -name "*.nwk"

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

In previous steps, we have manually edited `~/Scripts/withncbi/doc/mitochondrion_OG.md` and
generated `genus_OG.tsv`.

*D* between target and outgroup should be around **0.05**.

```bash
mkdir -p ~/data/organelle/mito/OG
cd ~/data/organelle/mito/OG

time bash ../mitochondrion.OG.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `; do
    echo "echo \"====> Processing ${d} <====\""
    echo bash ${d}/1_pair.sh;
    echo bash ${d}/2_rawphylo.sh;
    echo bash ${d}/3_multi.sh;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

## Self alignments

```bash
mkdir -p ~/data/organelle/mito/self
cd ~/data/organelle/mito/self

time bash ../mitochondrion.self.cmd.txt 2>&1 | tee log_cmd.txt

# Don't need 6_feature_cmd.sh 7_pair_stat.sh
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

# clean
find . -mindepth 1 -maxdepth 2 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | parallel -r rm -fr

```

