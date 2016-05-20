# Processing bacterial genomes species by species

The following command lines are about how I processed the plastid genomes of green plants.
Many tools of `taxon/` are used here, which makes a good example for users.

## Work flow.

```text
id ---> lineage ---> filtering ---> naming ---> strain_info.pl   ---> strain_bz.pl
                                      |                                 ^
                                      |-------> batch_get_seq.pl -------|
```

I'm sure there are no commas in names. So for convenient, don't use Text::CSV_XS.

## Init genome report database.

`db/README.md`

```bash
cd ~/Scripts/withncbi/db
perl gr_strains.pl -o prok_strains.csv
perl gr_db.pl --db gr_prok --file prok_strains.csv
```

Find valid species.

Got **160** species.

```bash
mkdir -p ~/data/bacteria/bac_summary
cd ~/data/bacteria/bac_summary

perl ~/Scripts/alignDB/util/query_sql.pl --db gr_prok -q '
        SELECT  species_id `#species_id`,
                COUNT(*) count,
                MAX(CHAR_LENGTH(code)) species_code
        FROM gr
        WHERE   1 = 1
        AND status LIKE "%Complete%"                # complete genomes
        AND species NOT LIKE "%Candidatus%"         # uncertainty classification
        AND taxonomy_id != species_id               # no strain ID
        AND organism_name NOT LIKE "%,%"            # avoid commas in names
        AND (chr IS NOT NULL OR LENGTH(CHR) > 0)    # has chromosome accession
        AND (LENGTH(wgs) = 0 OR wgs IS NULL)        # avoid bad assembly
        AND species_member > 2
        AND genus IS NOT NULL
        GROUP BY species_id
        HAVING count > 2 AND species_code > 0
        ORDER BY subgroup, species_id
    ' -o stdout \
    | cut -d ',' -f 1 \
    | grep -v "^#" \
    > bac.SPECIES.csv
cat bac.SPECIES.csv | wc -l
```

Expand species to strains. (Nested single quotes in bash should be '\'')

Got **1488** strains and **1604** chromosomes.


```bash
cd ~/data/bacteria/bac_summary

#echo '#strain_taxonomy_id,strain,species,genus,subgroup' > bac.STRAIN.csv

rm bac_strain.csv
cat bac.SPECIES.csv \
    | parallel --keep-order --no-run-if-empty -j 8 '
        perl ~/Scripts/alignDB/util/query_sql.pl --db gr_prok -q '\''
            SELECT  taxonomy_id `#strain_taxonomy_id`,
                    organism_name `strain`,
                    species,
                    genus,
                    subgroup,
                    chr
            FROM gr
            WHERE 1 = 1
            AND status LIKE "%Complete%"                # complete genomes
            AND species NOT LIKE "%Candidatus%"         # uncertainty classification
            AND taxonomy_id != species_id               # no strain ID
            AND organism_name NOT LIKE "%,%"            # avoid commas in names
            AND (chr IS NOT NULL OR LENGTH(CHR) > 0)    # has chromosome accession
            AND (LENGTH(wgs) = 0 OR wgs IS NULL)        # avoid bad assembly
            AND species_id = {}
        '\'' -o stdout
    ' \
    | grep -v "^#" \
    >> bac.STRAIN.csv
cat bac.STRAIN.csv | wc -l

cat bac.STRAIN.csv \
    | perl -nla -F"," -e '
        my @fields = splice @F, 0, 5;
        while (my $acc = shift @F) {
            $acc =~ s/"//g;
            $acc =~ s/\.\d+//;
            print join q{,}, (@fields, $acc);
        }
    ' \
    > bac.ACC.csv
cat bac.ACC.csv | wc -l
```

Create abbreviations.

```bash
cd ~/data/bacteria/bac_summary

echo '#strain_taxonomy_id,strain,species,genus,subgroup,accession,abbr' > bac.ABBR.csv
cat bac.ACC.csv \
    | grep -v '^#' \
    | perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "2,3,4" -s "," -m 0 \
    | sort -t',' -k5,5 -k4,4 -k3,3 -k7,7 \
    >> bac.ABBR.csv
```

## Download sequences and regenerate lineage information.

We don't rename sequences here, so the file has three columns.

And create `bac_ncbi.csv` with abbr names as taxon file.

```bash
mkdir -p ~/data/bacteria/bac_genomes
cd ~/data/bacteria/bac_genomes

echo "#strain_name,accession,strain_taxon_id" > bac_name_acc_id.csv
cat ../bac_summary/bac.ABBR.csv \
    | grep -v '^#' \
    | perl -nla -F"," -e 'print qq{$F[6],$F[5],$F[0]}' \
    | sort \
    >> bac_name_acc_id.csv

cat ../bac_summary/bac.ABBR.csv \
    | grep -v '^#' \
    | perl -nla -F"," -e 'print qq{$F[0],$F[6]}' \
    | uniq \
    | perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file bac_ncbi.csv

perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -f bac_name_acc_id.csv \
    -l ~/data/NCBI/genomes/Bacteria \
    -p 2>&1 \
    | tee bac_seq.log

# count downloaded sequences
find . -name "*.fasta" | wc -l
```

## Create alignment plans

We got **476** accessions.

Numbers for higher ranks are: 49 orders, 61 families, 119 genera and 467 species.

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

# 476
wc -l plastid.GENUS.csv

#   count every ranks
#   49 order.list.tmp
#   61 family.list.tmp
#  119 genus.list.tmp
#  467 species.list.tmp
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
    'BEGIN{ @ls = grep {/\S/} grep {!/^#/} path(q{~/Scripts/withncbi/doc/plastid_OG.md})->lines( { chomp => 1}); for (@ls) {@fs = split(/,/); $h{$fs[0]}= $fs[1];}  } if (exists $h{$F[0]}) { printf qq{%s\t%s\t%s\t%s\n}, $F[0] . q{_OG}, $F[1], $F[2], $h{$F[0]}; }' \
    > genus_OG.tsv

# every genera
echo -e "mkdir -p ~/data/organelle/plastid.working\ncd ~/data/organelle/plastid.working\n" > ../plastid.cmd.txt
cat genus.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl --seq_dir ~/data/organelle/plastid_genomes --taxon_file ~/data/organelle/plastid_genomes/plastid_ncbi.csv --parallel 8 \
    >> ../plastid.cmd.txt

# this is for finding outgroups
echo -e "mkdir -p ~/data/organelle/plastid_families\ncd ~/data/organelle/plastid_families\n" > ../plastid_families.cmd.txt
cat family.tsv \
    | perl -n -e '/,\w+,/ and print' \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl --seq_dir ~/data/organelle/plastid_genomes --taxon_file ~/data/organelle/plastid_genomes/plastid_ncbi.csv --parallel 8 \
    >> ../plastid_families.cmd.txt

# genera with outgroups
echo -e "mkdir -p ~/data/organelle/plastid_OG\ncd ~/data/organelle/plastid_OG\n" > ../plastid_OG.cmd.txt
cat genus_OG.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl --seq_dir ~/data/organelle/plastid_genomes --taxon_file ~/data/organelle/plastid_genomes/plastid_ncbi.csv --parallel 8 \
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
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
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

#----------------------------#
# Approach 2: step by step
#----------------------------#
# real_chr
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 1_real_chr.sh | sort `;do
    echo sh $f ;
    echo ;
done  > run_1.sh

# RepeatMasker
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 2_file_rm.sh | sort `;do
    echo sh $f ;
    echo ;
done  > run_2.sh

# pair
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 3_pair_cmd.sh | sort `;do
    echo sh $f ;
    echo ;
done  > run_3.sh

# rawphylo
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 4_rawphylo.sh | sort `;do
    echo sh $f ;
    echo ;
done  > run_4.sh

# multi cmd
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 5_multi_cmd.sh | sort `;do
    echo sh $f ;
    echo ;
done  > run_5.sh

# multi db
for f in `find . -mindepth 1 -maxdepth 2 -type f -name 7_multi_db_only.sh | sort `;do
    echo sh $f ;
    echo ;
done  > run_7.sh

cat run_1.sh | grep . | parallel -j 8 2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel -j 6 2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel -j 6 2>&1 | tee log_3.txt
cat run_4.sh | grep . | parallel -j 2 2>&1 | tee log_4.txt
cat run_5.sh | grep . | parallel -j 2 2>&1 | tee log_5.txt
cat run_7.sh | grep . | parallel -j 4 2>&1 | tee log_7.txt

#----------------------------#
# Charting on Windows
#----------------------------#
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    export d_base=`basename $d` ;
    echo "perl d:/Scripts/fig_table/collect_common_basic.pl    -d $d_base" ;
    echo "perl d:/Scripts/alignDB/stat/common_chart_factory.pl --replace diversity=divergence -i $d_base/$d_base.common.xlsx" ;
    echo "perl d:/Scripts/alignDB/stat/multi_chart_factory.pl  --replace diversity=divergence -i $d_base/$d_base.multi.xlsx" ;
    echo "perl d:/Scripts/alignDB/stat/gc_chart_factory.pl     --replace diversity=divergence -i $d_base/$d_base.gc.xlsx" ;
    echo ;
done  > run_chart.bat
perl -pi -e 's/\n/\r\n/g' run_chart.bat

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel --no-run-if-empty rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel --no-run-if-empty rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel --no-run-if-empty rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel --no-run-if-empty rm
```

### Self alignments.

```bash
cd ~/data/organelle/

perl -p -e 's/plastid\.working/plastid_self.working/g; s/strain_bz/strain_bz_self/g; s/(\-\-use_name)/\1 --length 1000 /g;' plastid.cmd.txt > plastid_self.cmd.txt

mkdir -p ~/data/organelle/plastid_self.working
cd ~/data/organelle/plastid_self.working

time sh ../plastid_self.cmd.txt 2>&1 | tee log_cmd.txt

# Don't need 6_feature_cmd.sh
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing $d <====\""
    echo sh $d/1_real_chr.sh ;
    echo sh $d/2_file_rm.sh ;
    echo sh $d/3_self_cmd.sh ;
    echo sh $d/4_proc_cmd.sh ;
    echo sh $d/5_circos_cmd ;
    echo sh $d/7_pair_stat.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

#----------------------------#
# Charting on Windows
#----------------------------#
for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    export d_base=`basename $d` ;
    echo "perl d:/Scripts/fig_table/collect_common_basic.pl    -d $d_base" ;
    echo "perl d:/Scripts/alignDB/stat/common_chart_factory.pl --replace diversity=divergence -i $d_base/${d_base}_paralog.common.xlsx" ;
    echo "perl d:/Scripts/alignDB/stat/gc_chart_factory.pl     --replace diversity=divergence -i $d_base/${d_base}_paralog.gc.xlsx" ;
    echo ;
done  > run_chart.bat
perl -pi -e 's/\n/\r\n/g' run_chart.bat

# clean
find . -mindepth 1 -maxdepth 2 -type d -name "*_raw" | parallel --no-run-if-empty rm -fr
find . -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | parallel --no-run-if-empty rm -fr

# clean mysql
find  /usr/local/var/mysql -type d -name "[A-Z]*" | parallel --no-run-if-empty rm -fr
```

### Alignments of families for outgroups.

The old prepare_run.sh

```bash
mkdir -p ~/data/organelle/plastid_families
cd ~/data/organelle/plastid_families

time sh ../plastid_families.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    echo "echo \"====> Processing $d <====\""
    echo sh $d/1_real_chr.sh ;
    echo sh $d/2_file_rm.sh ;
    echo sh $d/3_pair_cmd.sh ;
    echo sh $d/4_rawphylo.sh ;
    echo sh $d/5_multi_cmd.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    export d_base=`basename $d` ;
    echo "perl d:/Scripts/fig_table/collect_common_basic.pl    -d $d_base" ;
    echo ;
done  > run_chart.bat
perl -pi -e 's/\n/\r\n/g' run_chart.bat

# for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
#     export d_base=`basename $d` ;
#     export f_base="${d_base}/${d_base}_phylo/${d_base}" ;
#     if [ -f $f_base.nwk ]
#     then
#         echo $f_base ;  
#         nw_display -s -b 'visibility:hidden' $f_base.nwk > $f_base.svg ;
#     fi
# done

find . -type f -path "*_phylo*" -name "*.nwk"
```

Manually editing `~/Scripts/withncbi/doc/plastid_OG.md` and generate `genus_OG.tsv`.

*D* of outgroups should be around 0.05.

```bash
mkdir -p ~/data/organelle/plastid_OG
cd ~/data/organelle/plastid_OG

time sh ../plastid_OG.cmd.txt 2>&1 | tee log_cmd.txt

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

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `;do
    export d_base=`basename $d` ;
    echo "perl d:/Scripts/fig_table/collect_common_basic.pl    -d $d_base" ;
    echo "perl d:/Scripts/alignDB/stat/common_chart_factory.pl --replace diversity=divergence -i $d_base/$d_base.common.xlsx" ;
    echo "perl d:/Scripts/alignDB/stat/multi_chart_factory.pl  --replace diversity=divergence -i $d_base/$d_base.multi.xlsx" ;
    echo "perl d:/Scripts/alignDB/stat/gc_chart_factory.pl     --replace diversity=divergence -i $d_base/$d_base.gc.xlsx" ;
    echo ;
done  > run_chart.bat
perl -pi -e 's/\n/\r\n/g' run_chart.bat

find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel --no-run-if-empty rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel --no-run-if-empty rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel --no-run-if-empty rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel --no-run-if-empty rm

```

SNP t-test

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
| ------:  | :-----              | ------:      |
| 1163     | Anabaena            | 3            |
| 35823    | Arthrospira         | 2            |
| 1186     | Calothrix           | 3            |
| 102234   | Cyanobacterium      | 2            |
| 43988    | Cyanothece          | 6            |
| 33071    | Gloeobacter         | 2            |
| 1125     | Microcystis         | 2            |
| 1177     | Nostoc              | 4            |
| 1158     | Oscillatoria        | 2            |
| 1218     | Prochlorococcus     | 14           |
| 1129     | Synechococcus       | 20           |
| 1142     | Synechocystis       | 5            |
| 146785   | Thermosynechococcus | 2            |

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
| ------:    | :-----                     | ------:      |
| 118562     | Arthrospira platensis      | 2            |
| 1126       | Microcystis aeruginosa     | 2            |
| 1219       | Prochlorococcus marinus    | 12           |
| 32046      | Synechococcus elongatus    | 2            |
| 1148       | Synechocystis sp. PCC 6803 | 4            |

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
    echo sh $d/1_real_chr.sh ; \
    echo sh $d/2_file_rm.sh ; \
    echo sh $d/3_pair_cmd.sh ; \
    echo sh $d/4_rawphylo.sh ; \
    echo sh $d/5_multi_cmd.sh ; \
    echo sh $d/7_multi_db_only.sh ; \
    echo ; \
done  > runall.sh


find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel --no-run-if-empty rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel --no-run-if-empty rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel --no-run-if-empty rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel --no-run-if-empty rm
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

cat <<EOF > cmd_common_chart.tt
[% FOREACH item IN data -%]
[% chart = 'common' -%]
REM [% item.name %]
if not exist [% item.name %].[% chart %].xlsx goto skip[% item.name %]
perl d:/Scripts/alignDB/stat/[% chart %]_chart_factory.pl --replace diversity=divergence -i [% item.name %].[% chart %].xlsx
if not exist [% item.name %]_OG.[% chart %].xlsx goto skip[% item.name %]
perl d:/Scripts/alignDB/stat/[% chart %]_chart_factory.pl --replace diversity=divergence -i [% item.name %]_OG.[% chart %].xlsx
:skip[% item.name %]

[% END -%]
EOF

cat ~/data/organelle/plastid_summary/genus.tsv \
    | cut -f 1 \
    | TT_FILE=cmd_common_chart.tt perl -MTemplate -nl -e 'push @data, { name => $_, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, }) or die Template->error}' \
    > cmd_common_chart.bat

# Undre Windows
cd /d D:/data/organelle/plastid_summary/xlsx
cmd_common_chart.bat
```

### Genome list

Create `plastid.list.csv` from `plastid.GENUS.csv` with sequence lengths.

```bash
mkdir -p ~/data/organelle/plastid_summary/table
cd ~/data/organelle/plastid_summary/table

# manually set orders in `plastid_OG.md`
echo "genus" > genus_all.lst
perl -l -MPath::Tiny -e \
    'BEGIN{ @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_} map {s/,\w+//; $_} map {s/^###\s*//; $_} path(q{~/Scripts/withncbi/doc/plastid_OG.md})->lines( { chomp => 1}); } $fh; for (@ls) { (/^\s*$/ or /^##\s+/ or /^#\s+(\w+)/) and next; print $_}' \
    >> genus_all.lst

echo "genus,strain_abbr,accession,length" > length.tmp
find ~/data/organelle/plastid.working -type f -name "chr.sizes" \
    | xargs perl -nl -e \
    '$_ =~ s/\t/\,/; $ARGV =~ /working\/(\w+)\/(\w+)\//; print qq{$1,$2,$_}' \
    >> length.tmp

cat length.tmp \
    | perl -nl -a -F',' -MPath::Tiny -e \
    'BEGIN{ @ls = path(q{genus_all.lst})->lines( { chomp => 1}); $o{$ls[$_]} = $_ for (0 .. $#ls); } print qq{$_,$o{$F[0]}};' \
    | sort -n -t, -k5,5 \
    | cut -d',' -f 1-4 \
    > length_ordered.tmp

perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t length_ordered.tmp -m ~/data/organelle/plastid_summary/plastid.GENUS.csv -f 2 -f2 1 --concat --stdout \
    | perl -nl -a -F"," -e 'print qq{$F[12],$F[10],$F[8],$F[6],$F[4],$F[5],$F[3]}' \
    >  plastid.list.csv

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
cat <<EOF > Table_alignment.tt
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
    | grep -v "^genus" \
    | TT_FILE=Table_alignment.tt perl -MTemplate -nl -e 'push @data, { name => $_, file => qq{$_.common.xlsx}, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, }) or die Template->error}' \
    > Table_alignment_all.yml

# Under Windows
cd /d D:/data/organelle/plastid_summary/xlsx
perl d:/Scripts/fig_table/excel_table.pl -i Table_alignment_all.yml
perl d:/Scripts/fig_table/xlsx2xls.pl -d Table_alignment_all.xlsx --csv

# Back to Mac
cd ~/data/organelle/plastid_summary/xlsx
cp -f Table_alignment_all.xlsx ~/data/organelle/plastid_summary/table
perl -pi -e 's/\r\n/\n/g' Table_alignment_all.csv
cp -f Table_alignment_all.csv ~/data/organelle/plastid_summary/table

cd ~/data/organelle/plastid_summary/table

echo "Genus,avg_size" > genus_avg_size.csv
cat plastid.list.csv \
    | grep -v "#" \
    | perl -nl -a -F"," -e \
    '$count{$F[2]}++; $sum{$F[2]} += $F[6]; END {for $k (sort keys %count) { printf qq{%s,%d\n}, $k, $sum{$k}/$count{$k}; } }' \
    >> genus_avg_size.csv

perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t Table_alignment_all.csv -m genus_avg_size.csv -f 0 -f2 0 --concat --stdout \
    > Table_alignment_all.1.csv

echo "Genus,coverage" > genus_coverage.csv
cat Table_alignment_all.1.csv \
    | perl -nl -a -F',' -e '$F[7] =~ /\D+/ and next; $c = $F[2] * 1000 * 1000 / $F[7]; print qq{$F[0],$c};' \
    >> genus_coverage.csv

perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t Table_alignment_all.1.csv -m genus_coverage.csv -f 0 -f2 0 --concat --stdout \
    > Table_alignment_all.2.csv

echo "Genus,indels" > genus_indels.csv
cat Table_alignment_all.2.csv \
    | perl -nl -a -F',' -e '$F[7] =~ /\D+/ and next; $c = $F[3] / 100 * $F[2] * 1000 * 1000; print qq{$F[0],$c};' \
    >> genus_indels.csv

perl ~/Scripts/alignDB/util/merge_csv.pl \
    -t Table_alignment_all.2.csv -m genus_indels.csv -f 0 -f2 0 --concat --stdout \
    > Table_alignment_all.3.csv

cat Table_alignment_all.3.csv \
    | cut -d',' -f 1-6,8,10,12 \
    > Table_alignment_for_filter.csv

# real filter
cat Table_alignment_for_filter.csv \
    | perl -nl -a -F',' -e '$F[6] =~ /\D+/ and next; print $F[0] if ($F[7] < 0.4 or $F[8] < 100 or $F[4] > 0.2);' \
    > genus_exclude.lst

grep -v -Fx -f genus_exclude.lst genus_all.lst > genus.lst

rm ~/data/organelle/plastid_summary/table/Table_alignment_all.[0-9].csv
rm ~/data/organelle/plastid_summary/table/genus*csv

#
cd ~/data/organelle/plastid_summary/xlsx
cat ~/data/organelle/plastid_summary/table/genus.lst \
    | grep -v "^genus" \
    | TT_FILE=Table_alignment.tt perl -MTemplate -nl -e 'push @data, { name => $_, file => qq{$_.common.xlsx}, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, }) or die Template->error}' \
    > Table_alignment.yml

# Under Windows
cd /d D:/data/organelle/plastid_summary/xlsx
perl d:/Scripts/fig_table/excel_table.pl -i Table_alignment.yml
perl d:/Scripts/fig_table/xlsx2xls.pl -d Table_alignment.xlsx --csv

# Mac
cp -f ~/data/organelle/plastid_summary/xlsx/Table_alignment.xlsx ~/data/organelle/plastid_summary/table
perl -pi -e 's/\r\n/\n/g' ~/data/organelle/plastid_summary/xlsx/Table_alignment.csv
cp -f ~/data/organelle/plastid_summary/xlsx/Table_alignment.csv ~/data/organelle/plastid_summary/table
```

### Groups

```bash
mkdir -p ~/data/organelle/plastid_summary/group
cd ~/data/organelle/plastid_summary/group

perl -l -MPath::Tiny -e \
'BEGIN{ @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_} map {s/,\w+//; $_} map {s/^###\s*//; $_} path(q{~/Scripts/withncbi/doc/plastid_OG.md})->lines( { chomp => 1}); } $fh; for (@ls) { (/^\s*$/ or /^##\s+/) and next; if (/^#\s+(\w+)/) {$fh = path("$1.txt")->openw; next;} else {print {$fh} $_}}'

grep -Fx -f ~/data/organelle/plastid_summary/table/genus.lst Angiosperm.txt > Angiosperm.lst
grep -Fx -f ~/data/organelle/plastid_summary/table/genus.lst Gymnosperm.txt > Gymnosperm.lst

find . -type f -name "*.txt" \
    | xargs cat \
    | grep -Fx -f ~/data/organelle/plastid_summary/table/genus.lst \
    > Others.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nl -a -F',' -e '$F[1] > 0.05 and print $F[0];' \
    > group_4.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nl -a -F',' -e '$F[1] > 0.02 and $F[1] <= 0.05 and print $F[0];' \
    > group_3.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nl -a -F',' -e '$F[1] > 0.005 and $F[1] <= 0.02 and print $F[0];' \
    > group_2.lst

cat ~/data/organelle/plastid_summary/table/Table_alignment.csv \
    | cut -d, -f 1,5 \
    | perl -nl -a -F',' -e '$F[1] <= 0.005 and print $F[0];' \
    > group_1.lst

rm *.txt
```

NCBI Taxonomy tree

```bash
cd ~/data/organelle/plastid_summary/group

cat ~/data/organelle/plastid_summary/table/genus.lst \
    | perl -e '@ls = <>; $str = qq{bp_taxonomy2tree.pl \\\n}; for (@ls) {chomp;$str .= qq{-s $_ \\\n}}  $str .= qq{-e \n}; print $str' \
    > genera_tree.sh
sh genera_tree.sh > genera.tree

```

### Phylogenic trees of each genus with outgroup

```bash
mkdir -p ~/data/organelle/plastid_summary/trees

cat ~/Scripts/withncbi/doc/plastid_OG.md \
    | grep -v "^#" | grep . \
    | cut -d',' -f 1 \
    > ~/data/organelle/plastid_summary/trees/list.txt

find ~/data/organelle/plastid_OG -type f -path "*_phylo*" -name "*.nwk" \
    | parallel -j 1 cp {} ~/data/organelle/plastid_summary/trees
```

### d1, d2

`collect_excel.pl`

```bash
cd ~/data/organelle/plastid_summary/xlsx

cat <<EOF > cmd_collect_d1_d2.tt
perl d:/Scripts/fig_table/collect_excel.pl [% FOREACH item IN data -%] -f [% item.name %].common.xlsx -s d1_pi_gc_cv -n [% item.name %] [% END -%] -o cmd_plastid_d1.xlsx

perl d:/Scripts/fig_table/collect_excel.pl [% FOREACH item IN data -%] -f [% item.name %].common.xlsx -s d2_pi_gc_cv -n [% item.name %] [% END -%] -o cmd_plastid_d2.xlsx

perl d:/Scripts/fig_table/collect_excel.pl [% FOREACH item IN data -%] -f [% item.name %].common.xlsx -s d1_comb_pi_gc_cv -n [% item.name %] [% END -%] -o cmd_plastid_d1_comb.xlsx

perl d:/Scripts/fig_table/collect_excel.pl [% FOREACH item IN data -%] -f [% item.name %].common.xlsx -s d2_comb_pi_gc_cv -n [% item.name %] [% END -%] -o cmd_plastid_d2_comb.xlsx

EOF

cat ~/data/organelle/plastid_summary/table/genus.lst \
    | TT_FILE=cmd_collect_d1_d2.tt perl -MTemplate -nl -e 'push @data, { name => $_, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, }) or die Template->error}' \
    > cmd_collect_d1_d2.bat

# Undre Windows
cd /d D:/data/organelle/plastid_summary/xlsx
cmd_collect_d1_d2.bat
```

`ofg_chart.pl`

```bash
mkdir -p ~/data/organelle/plastid_summary/fig

cd ~/data/organelle/plastid_summary/xlsx

cat <<EOF > cmd_chart_d1_d2.tt
perl d:/Scripts/fig_table/ofg_chart.pl -i cmd_plastid_d1.xlsx -xl "Distance to indels ({italic(d)[1]})" -yl "Nucleotide divergence ({italic(D)})" -xr "A2:A8" -yr "B2:B8"  --y_min 0.0 --y_max [% y_max %] -x_min 0 -x_max 5 -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" -rs "NON_EXIST" --postfix [% postfix %] --style_dot -ms

perl d:/Scripts/fig_table/ofg_chart.pl -i cmd_plastid_d1_comb.xlsx -xl "Distance to indels ({italic(d)[1]})" -yl "Nucleotide divergence ({italic(D)})" -xr "A2:A8" -yr "B2:B8"  --y_min 0.0 --y_max [% y_max %] -x_min 0 -x_max 5 -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" -rs "NON_EXIST" --postfix [% postfix %] --style_dot

perl d:/Scripts/fig_table/ofg_chart.pl -i cmd_plastid_d2.xlsx -xl "Reciprocal of indel density ({italic(d)[2]})" -yl "Nucleotide divergence ({italic(D)})" -xr "A2:A23" -yr "B2:B23"  --y_min 0.0 --y_max [% y_max2 %] -x_min 0 -x_max 20 -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" -rs "NON_EXIST" --postfix [% postfix %] --style_dot  -ms

perl d:/Scripts/fig_table/ofg_chart.pl -i cmd_plastid_d2_comb.xlsx -xl "Reciprocal of indel density ({italic(d)[2]})" -yl "Nucleotide divergence ({italic(D)})" -xr "A2:A23" -yr "B2:B23"  --y_min 0.0 --y_max [% y_max2 %] -x_min 0 -x_max 20 -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" -rs "NON_EXIST" --postfix [% postfix %] --style_dot

EOF

cat ~/data/organelle/plastid_summary/group/group_1.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e 'push @data, { name => $_, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, y_max => 0.01, y_max2 => 0.01, postfix => q{group_1}, }) or die Template->error}' \
    > cmd_chart_group_1.bat

cat ~/data/organelle/plastid_summary/group/group_2.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e 'push @data, { name => $_, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, y_max => 0.03, y_max2 => 0.04, postfix => q{group_2}, }) or die Template->error}' \
    > cmd_chart_group_2.bat

cat ~/data/organelle/plastid_summary/group/group_3.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e 'push @data, { name => $_, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, y_max => 0.05, y_max2 => 0.05, postfix => q{group_3}, }) or die Template->error}' \
    > cmd_chart_group_3.bat

cat ~/data/organelle/plastid_summary/group/group_4.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e 'push @data, { name => $_, }; END{$tt = Template->new; $tt->process($ENV{TT_FILE}, { data => \@data, y_max => 0.15, y_max2 => 0.15, postfix => q{group_4}, }) or die Template->error}' \
    > cmd_chart_group_4.bat

# Undre Windows
cd /d D:/data/organelle/plastid_summary/xlsx
cmd_chart_group_1.bat
cmd_chart_group_2.bat
cmd_chart_group_3.bat
cmd_chart_group_4.bat

# Mac
rm ~/data/organelle/plastid_summary/xlsx/*.csv
cp ~/data/organelle/plastid_summary/xlsx/*.pdf ~/data/organelle/plastid_summary/fig

# Coreldraw doesn't play well with computer modern fonts (latex math).
# perl ~/Scripts/fig_table/tikz_chart.pl -i cmd_plastid_d1_A2A8_B2B8.group_1.csv -xl 'Distance to indels ($d_1$)' -yl 'Nucleotide divergence ($D$)' --y_min 0.0 --y_max 0.01 -x_min 0 -x_max 5 --style_dot --pdf

```
