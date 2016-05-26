# Processing bacterial genomes species by species

The following command lines are about how I processed the plastid genomes of green plants.
Many tools of `taxon/` are used here, which makes a good example for users.

## Init genome report database.

`db/README.md`

```bash
cd ~/Scripts/withncbi/db
perl gr_strains.pl -o prok_strains.csv
perl gr_db.pl --db gr_prok --file prok_strains.csv
rm prok_strains.csv
```

Find valid species.

Got **159** species.

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
        HAVING count > 2 AND species_code > 0       # having enough and representative member
        ORDER BY subgroup, species_id
    ' -o stdout \
    | cut -d ',' -f 1 \
    | grep -v "^#" \
    > bac.SPECIES.csv
cat bac.SPECIES.csv | wc -l
```

Expand species to strains. (Nested single quotes in bash should be '\'')

Got **1470** strains.


```bash
cd ~/data/bacteria/bac_summary

rm bac.STRAIN.csv
cat bac.SPECIES.csv \
    | parallel --keep-order --no-run-if-empty -j 8 '
        perl ~/Scripts/alignDB/util/query_sql.pl --db gr_prok -q '\''
            SELECT  taxonomy_id `#strain_taxonomy_id`,
                    organism_name `strain`,
                    species,
                    genus,
                    subgroup,
                    `code`,
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
```

Create abbreviations.

```bash
cd ~/data/bacteria/bac_summary

echo '#strain_taxonomy_id,strain,species,genus,subgroup,code,accession,abbr' > bac.ABBR.csv
cat bac.STRAIN.csv \
    | grep -v '^#' \
    | perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "2,3,4" -s "," -m 0 --shortsub \
    | sort -t',' -k5,5 -k4,4 -k3,3 -k6,6 \
    >> bac.ABBR.csv
```

## Download sequences and regenerate lineage information.

We don't rename sequences here, so the file has three columns. **1587** accessions.

And create `bac_ncbi.csv` with abbr names as taxon file.

```bash
mkdir -p ~/data/bacteria/bac_genomes
cd ~/data/bacteria/bac_genomes

cat ../bac_summary/bac.ABBR.csv \
    | grep -v '^#' \
    | perl -nla -F"," -e 'print qq{$F[0],$F[7]}' \
    | uniq \
    | perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file bac_ncbi.csv

echo "#strain_name,accession,strain_taxon_id" > bac_name_acc_id.csv
cat ../bac_summary/bac.ABBR.csv \
    | grep -v '^#' \
    | perl -nla -F"," -e '
        my $acc = $F[6];
        $acc =~ s/"//g;
        $acc =~ s/\.\d+//g;
        for my $s (split /\|/, $acc) {
            print qq{$F[7],$s,$F[0]};
        }
    ' \
    | sort \
    >> bac_name_acc_id.csv
cat bac_name_acc_id.csv | wc -l

perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -f bac_name_acc_id.csv \
    -l ~/data/NCBI/genomes/Bacteria \
    -p 2>&1 \
    | tee bac_seq.log

# count downloaded sequences
find . -name "*.fasta" | wc -l
```

## Create alignment plans

Numbers for higher ranks are: 16 subgroups, 80 genera and 159 species.

```bash
cd ~/data/bacteria/bac_summary/

# count every ranks
#  16 subgroup.list.tmp
#  80 genus.list.tmp
# 159 species.list.tmp
cat bac.ABBR.csv | grep -v '^#'| cut -d',' -f 3 | sort | uniq > species.list.tmp
cat bac.ABBR.csv | grep -v '^#'| cut -d',' -f 4 | sort | uniq > genus.list.tmp
cat bac.ABBR.csv | grep -v '^#'| cut -d',' -f 5 | sort | uniq > subgroup.list.tmp
wc -l subgroup.list.tmp genus.list.tmp species.list.tmp

rm *.tmp
```

Create `bac_target_OG.md` for picking target and outgroup.

Manually edit it then move to `~/Scripts/withncbi/doc/bac_target_OG.md`.

```bash
cd ~/data/bacteria/bac_summary

cat bac.ABBR.csv \
    | grep -v "^#" \
    | perl -na -F"," -e '
        BEGIN{
            ($subgroup, $genus, $species,) = (q{}, q{}, q{});
        }

        next if ! $F[5]; # code
        chomp for @F;

        if ($F[4] ne $subgroup) {
            $subgroup = $F[4];
            printf qq{\n# %s\n}, $subgroup;
        }
        if ($F[3] ne $genus) {
            $genus = $F[3];
            printf qq{## %s\n}, $genus;
        }
        $F[2] =~ s/\W+/_/g;
        if ($F[2] ne $species) {
            $species = $F[2];
        }
        printf qq{%s,%s,%s\n}, $species, $F[7], $F[5];
    ' \
    > bac_target_OG.md
```

Create alignments without outgroups.

```text
ABBR.csv
#strain_taxonomy_id,strain,species,genus,subgroup,code,accession,abbr
```

```bash
cd ~/data/bacteria/bac_summary

# 159, same as previous count
cat ~/Scripts/withncbi/doc/bac_target_OG.md \
    | grep -v '^#' \
    | grep -E '\S+' \
    | wc -l

# tab-separated
# name  t   qs
cat bac.ABBR.csv \
    | perl -nl -a -F"," -MPath::Tiny -e '
        BEGIN{
            $name = q{};
            %id_of = ();
            %h = ();
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{~/Scripts/withncbi/doc/bac_target_OG.md})->lines({chomp => 1});
            for (@ls) {
                @fs = split(/,/);
                $h{$fs[0]}= $fs[1];
            }
            undef @ls;
        }

        chomp for @F;
        $F[2] =~ s/\W+/_/g;
        if ($F[2] ne $name) {
            if ($name) {
                if (exists $h{$name}) {
                    my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
                    my $t = $h{$name};
                    my $qs = join(q{,}, grep {$_ ne $h{$name}} @s);
                    printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
                }
            }
            $name = $F[2];
            %id_of = ();
        }
        $id_of{$F[7]} = $F[0]; # same strain multiple chromosomes collapsed here

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = $h{$name};
            my $qs = join(q{,}, grep {$_ ne $h{$name}} @s);
            printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
        }' \
    > species.tsv

# every species
cd ~/data/bacteria/bac_summary

echo -e "mkdir -p ~/data/bacteria/bac.working \ncd ~/data/bacteria/bac.working\n" > ../bac.cmd.txt
cat species.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl \
        --seq_dir ~/data/bacteria/bac_genomes/ \
        --csv_taxon ~/data/bacteria/bac_genomes/bac_ncbi.csv \
        --parallel 8 \
    >> ../bac.cmd.txt

echo -e "mkdir -p ~/data/bacteria/bac.working \ncd ~/data/bacteria/bac.working\n" > ../bac.redo.cmd.txt
cat species.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl \
        --csv_taxon ~/data/bacteria/bac_genomes/bac_ncbi.csv \
        --parallel 8 \
    >> ../bac.redo.cmd.txt
```

Align all representative strains of every genera.

**39** genera.

```bash
cd ~/data/bacteria/bac_summary

cat bac.ABBR.csv \
    | grep -v "^#" \
    | perl -na -F"," -e '
        BEGIN{
            $name = q{};
            %id_of = ();
        }

        chomp for @F;
        next if ! $F[5]; # code
        my $genus = $F[3];
        $genus =~ s/\W+/_/g;
        if ($genus ne $name) {
            if ($name) {
                # sort by taxonomy_id
                my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
                my $t = shift @s;
                my $qs = join(q{,}, @s);
                printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
            }
            $name = $genus;
            %id_of = ();
        }
        $id_of{$F[7]} = $F[0]; # multiple chromosomes collapsed here

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = shift @s;
            my $qs = join(q{,}, @s);
            printf qq{%s\t%s\t%s\n}, $name, $t, $qs;
        }
    ' \
    | perl -nl -e '/\w+\t\w+\t\w+/ and print' \
    > genus.tsv
cat genus.tsv | wc -l

# every species
echo -e "mkdir -p ~/data/bacteria/bac.genus.working \ncd ~/data/bacteria/bac.genus.working" > ../bac.genus.cmd.txt
cat genus.tsv \
    | perl ~/Scripts/withncbi/taxon/cmd_template.pl \
        --seq_dir ~/data/bacteria/bac_genomes/ \
        --csv_taxon ~/data/bacteria/bac_genomes/bac_ncbi.csv \
        --parallel 8 \
    >> ../bac.genus.cmd.txt
```


## Aligning

### Batch running for groups

The old prepare_run.sh

```bash
mkdir -p ~/data/bacteria/bac.working
cd ~/data/bacteria/bac.working

bash ../bac.cmd.txt 2>&1 | tee log_cmd.txt
# bash ../bac.redo.cmd.txt 2>&1 | tee log_redo_cmd.txt # skip real_chr and repeatmasker

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
cat run_1.sh | grep . | parallel --no-run-if-empty -j 16 2>&1 | tee log_1.txt
cat run_2.sh | grep . | parallel --no-run-if-empty -j 8  2>&1 | tee log_2.txt
cat run_3.sh | grep . | parallel --no-run-if-empty -j 16 2>&1 | tee log_3.txt
cat run_4.sh | grep . | parallel --no-run-if-empty -j 4  2>&1 | tee log_4.txt
cat run_5.sh | grep . | parallel --no-run-if-empty -j 4  2>&1 | tee log_5.txt
cat run_7.sh | grep . | parallel --no-run-if-empty -j 16 2>&1 | tee log_7.txt

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel --no-run-if-empty rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel --no-run-if-empty rm -fr

find . -mindepth 1 -maxdepth 4 -type f -name "*.phy" | parallel --no-run-if-empty rm
find . -mindepth 1 -maxdepth 4 -type f -name "*.phy.reduced" | parallel --no-run-if-empty rm
```

### Alignments of genera for outgroups.

```bash
mkdir -p ~/data/bacteria/bac.genus.working
cd ~/data/bacteria/bac.genus.working

time bash ../bac.genus.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `; do
    echo "echo \"====> Processing ${d} <====\""
    echo sh ${d}/1_real_chr.sh ;
    echo sh ${d}/2_file_rm.sh ;
    echo sh ${d}/3_pair_cmd.sh ;
    echo sh ${d}/4_rawphylo.sh ;
    echo sh ${d}/5_multi_cmd.sh ;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

# Clean
```

## Summary

### Copy xlsx files

```bash
mkdir -p ~/data/bacteria/bac_summary/xlsx
cd ~/data/bacteria/bac_summary/xlsx

find  ~/data/bacteria/bac.working -type f -name "*.common.xlsx" \
    | grep -v "vs[A-Z]" \
    | parallel cp {} .

```

### Genome list

Create `bac.list.csv` from `bac.ABBR.csv` with sequence lengths.

```bash
mkdir -p ~/data/bacteria/bac_summary/table
cd ~/data/bacteria/bac_summary/table

# manually set orders in `bac_target_OG.md`
echo "#species" > species_all.lst
perl -l -MPath::Tiny -e '
    BEGIN {
        @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_} 
            map {s/,\w+//; $_} 
            map {s/^###\s*//; $_} 
            path(q{~/Scripts/withncbi/doc/bac_target_OG.md})->lines({chomp => 1});
    }
    
    for (@ls) { 
        (/^\s*$/ or /^##\s+/ or /^#\s+(\w+)/) and next; 
        print $_;
    }
    ' \
    >> species_all.lst

echo "#abbr,species,accession,length" > length.tmp
find ~/data/bacteria/bac.working -type f -name "chr.sizes" \
    | parallel --jobs 1 --keep-order --no-run-if-empty '
        perl -nl -e '\''
            BEGIN {
                my %l = ();
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

echo "#abbr,subgroup,genus,species,taxon_id" > abbr.tmp
cat ~/data/bacteria/bac_summary/bac.ABBR.csv \
    | grep -v "^#" \
    | perl -nla -F"," -e 'print qq{$F[7],$F[4],$F[3],$F[2],$F[0]}' \
    >> abbr.tmp

# #abbr,species,accession,length,subgroup,genus,species,taxon_id
cat length.tmp abbr.tmp \
    | perl ~/Scripts/withncbi/util/merge_csv.pl \
    -f 0 --concat -o stdout \
    | perl -nl -a -F"," -e 'print qq{$F[4],$F[5],$F[6],$F[0],$F[7],$F[2],$F[3]}' \
    > list.tmp

echo "#subgroup,genus,species,abbr,taxon_id,accession,length" > bac.list.csv
cat list.tmp \
    | grep -v "#" \
    | perl -nl -a -F',' -MPath::Tiny -e '
        BEGIN{ 
            @ls = path(q{species_all.lst})->lines({ chomp => 1}); 
            $o{$ls[$_]} = $_ for (0 .. $#ls); 
        } 
        $F[2] =~ s/\s+/_/g;
        print qq{$_,$o{$F[2]}};
    ' \
    | sort -n -t, -k8,8 \
    | cut -d',' -f 1-7 \
    >> bac.list.csv

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
