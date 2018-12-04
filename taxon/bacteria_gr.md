# Processing bacterial genomes species by species

[TOC]: # " "
- [Batch running for groups](#batch-running-for-groups)
- [Alignments of genera for outgroups.](#alignments-of-genera-for-outgroups)
- [Copy xlsx files](#copy-xlsx-files)
- [Genome list](#genome-list)
- [Genome alignment statistics](#genome-alignment-statistics)
- [sep_chart of d1, d2](#sep_chart-of-d1-d2)
- [CorelDRAW GC charts](#coreldraw-gc-charts)


# Init genome report database.

* Create database by following steps in
  [`db/README.md`](https://github.com/wang-q/withncbi/blob/master/db/README.md#genome-reports)

* Find valid species.

    * Got **207** species.
    * 40 species have no `species_code`

```bash
mkdir -p ~/data/bacteria/summary
cd ~/data/bacteria/summary

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
        HAVING count > 2 #AND species_code > 0       # having enough and representative member
        ORDER BY subgroup, species_id
    ' -o stdout |
    cut -d ',' -f 1 |
    grep -v "^#" \
    > bac.SPECIES.csv
cat bac.SPECIES.csv | wc -l

```

* Expand species to strains. (Nested single quotes in bash should be '\'')

    Got **2091** strains.

```bash
cd ~/data/bacteria/summary

cat bac.SPECIES.csv |
    parallel --keep-order -r -j 8 '
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
            ORDER BY released_date                      # oldest first
        '\'' -o stdout |
            grep -v "^#" |
            datamash transpose -t, |
            perl -MList::MoreUtils -e '\''
                my @lines = <>;
                if ($lines[5] =~ /^,+$/) {
                    $lines[5] = q{first} . $lines[5];
                }
                print $_ for @lines;
            '\'' |
            datamash transpose -t,
    ' |
    grep -v "^#" \
    > bac.STRAIN.csv
cat bac.STRAIN.csv | wc -l

```

* Create abbreviations.

```bash
cd ~/data/bacteria/summary

echo '#strain_taxonomy_id,strain,species,genus,subgroup,code,accession,abbr' > bac.ABBR.csv
cat bac.STRAIN.csv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "2,3,4" -s "," -m 0 --shortsub |
    sort -t',' -k5,5 -k4,4 -k3,3 -k6,6 \
    >> bac.ABBR.csv

```

# Download sequences and regenerate lineage information.

We don't rename sequences here, so the file has three columns. **2231** accessions.

And create `bac_ncbi.csv` with abbr names as taxon file.

```bash
mkdir -p ~/data/bacteria/GENOMES
cd ~/data/bacteria/GENOMES

cat ../summary/bac.ABBR.csv |
    grep -v '^#' |
    perl -nla -F"," -e 'print qq{$F[0],$F[7]}' |
    uniq |
    perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file bac_ncbi.csv

echo "#strain_name,accession,strain_taxon_id" > bac_name_acc_id.csv
cat ../summary/bac.ABBR.csv |
    grep -v '^#' |
    perl -nla -F"," -e '
        my $acc = $F[6];
        $acc =~ s/"//g;
        $acc =~ s/\.\d+//g;
        for my $s (split /\|/, $acc) {
            print qq{$F[7],$s,$F[0]};
        }
    ' |
    sort \
    >> bac_name_acc_id.csv
cat bac_name_acc_id.csv | wc -l

perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
    -f bac_name_acc_id.csv \
    -l ~/data/NCBI/genomes/Bacteria \
    2>&1 |
    tee bac_seq.log

# count downloaded sequences
find . -maxdepth 2 -name "*.fa" | wc -l
find . -maxdepth 1 -type d | wc -l

```

# Prepare sequences for lastz

```bash
cd ~/data/bacteria/GENOMES

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

Numbers for higher ranks are: 18 subgroups, 104 genera and 207 species.

```bash
cd ~/data/bacteria/summary/

# count every ranks
#  18 subgroup.list.tmp
# 104 genus.list.tmp
# 207 species.list.tmp
cat bac.ABBR.csv | grep -v '^#'| cut -d',' -f 3 | sort | uniq > species.list.tmp
cat bac.ABBR.csv | grep -v '^#'| cut -d',' -f 4 | sort | uniq > genus.list.tmp
cat bac.ABBR.csv | grep -v '^#'| cut -d',' -f 5 | sort | uniq > subgroup.list.tmp
wc -l subgroup.list.tmp genus.list.tmp species.list.tmp

rm *.tmp
```

Exclude diverged strains.

* 391904, Bifidobacterium longum subsp. infantis ATCC 15697 = JCM 1222 = DSM 20088, 2008-11-20,
  Complete Genome,
* 553190, Gardnerella vaginalis 409-05, 2010-01-07, Complete Genome, COM
* 1386087, Neisseria meningitidis LNP21362, 2015-01-07, Complete Genome,
* 935590, Neisseria meningitidis M0579, 2015-06-19, Complete Genome,
* 1415774, Clostridium botulinum 202F, 2014-12-04, Complete Genome,
* 508767, Clostridium botulinum E3 str. Alaska E43, 2008-05-16, Complete Genome,
* 929506, Clostridium botulinum BKT015925, 2011-04-18, Complete Genome,
* 935198, Clostridium botulinum B str. Eklund 17B (NRP), 2008-05-07, Complete Genome,
* 869303, Streptococcus pneumoniae SPN034156, 2010-07-29, Complete Genome,
* 869311, Streptococcus pneumoniae SPN032672, 2010-07-29, Complete Genome,
* 869312, Streptococcus pneumoniae SPN033038, 2010-07-29, Complete Genome,
* 261317, Buchnera aphidicola (Cinara tujafilina), 2011-06-09, Complete Genome,
* 372461, Buchnera aphidicola BCc, 2006-10-18, Complete Genome,

```sql
SELECT taxonomy_id, organism_name, released_date, status, code
FROM gr_prok.gr
WHERE species = "Gluconobacter oxydans" 
and status NOT IN ('Contig', 'Scaffold')
ORDER BY released_date, status, code

```

```bash
cd ~/data/bacteria/summary

cat bac.ABBR.csv |
    grep -v "391904," |
    grep -v "553190," |
    grep -v "1386087," |
    grep -v "935590," |
    grep -v "1415774," |
    grep -v "508767," |
    grep -v "929506," |
    grep -v "935198," |
    grep -v "869303," |
    grep -v "869311," |
    grep -v "869312," |
    grep -v "261317," |
    grep -v "372461," \
    > bac.WORKING.csv

```

Create `bac_target_OG.md` for picking target and outgroup.

Manually edit it then move to `~/Scripts/withncbi/doc/bac_target_OG.md`.

```bash
cd ~/data/bacteria/summary

cat bac.ABBR.csv |
    grep -v "^#" |
    perl -na -F"," -e '
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
bac.WORKING.csv
#strain_taxonomy_id,strain,species,genus,subgroup,code,accession,abbr
```

```bash
cd ~/data/bacteria/summary

# 207
cat ~/Scripts/withncbi/doc/bac_target_OG.md |
    grep -v '^#' |
    grep -E '\S+' |
    wc -l

# tab-separated
# name  t   qs
cat bac.WORKING.csv |
    perl -nl -a -F"," -MPath::Tiny -e '
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
cd ~/data/bacteria/summary

cat <<'EOF' > egaz_template_multi.tt

# [% name %]
egaz template \
    ~/data/bacteria/GENOMES/[% t %] \
[% FOREACH q IN qs -%]
    ~/data/bacteria/GENOMES/[% q %] \
[% END -%]
[% IF o -%]
    ~/data/bacteria/GENOMES/[% o %] \
    --outgroup [% o %] \
[% END -%]
    --multi -o [% name %] \
    --taxon ~/data/bacteria/GENOMES/bac_ncbi.csv \
    --rawphylo --parallel 8 -v

EOF

# every species
echo "mkdir -p ~/data/bacteria/species"  > ../bac.cmd.txt
echo "cd       ~/data/bacteria/species" >> ../bac.cmd.txt

cat species.tsv |
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
    >> ../bac.cmd.txt

cat <<'EOF' > egaz_templates_self.tt

# [% name %]
egaz template \
    ~/data/bacteria/GENOMES/[% t %] \
[% FOREACH q IN qs -%]
    ~/data/bacteria/GENOMES/[% q %] \
[% END -%]
    --self -o [% name %] \
    --taxon ~/data/bacteria/GENOMES/bac_ncbi.csv \
    --circos --aligndb --parallel 8 -v

EOF

# every species
echo "mkdir -p ~/data/bacteria/self"  > ../bac.self.cmd.txt
echo "cd       ~/data/bacteria/self" >> ../bac.self.cmd.txt

cat species.tsv |
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
    >> ../bac.self.cmd.txt

```

Align all representative strains of every genera.

**44** genera.

```bash
cd ~/data/bacteria/summary

cat bac.WORKING.csv |
    grep -v "^#" |
    perl -na -F"," -e '
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
    ' |
    perl -nl -e '/\w+\t\w+\t\w+/ and print' \
    > genus.tsv
cat genus.tsv | wc -l

# every genera
echo "mkdir -p ~/data/bacteria/genus"  > ../bac.genus.cmd.txt
echo "cd       ~/data/bacteria/genus" >> ../bac.genus.cmd.txt

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
    >> ../bac.genus.cmd.txt

```

# Aligning

## Batch running for groups

```bash
mkdir -p ~/data/bacteria/species
cd ~/data/bacteria/species

bash ../bac.cmd.txt 2>&1 | tee log_cmd.txt

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

#----------------------------#
# Clean
#----------------------------#
find . -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

## Alignments of genera for outgroups.

```bash
mkdir -p ~/data/bacteria/genus
cd ~/data/bacteria/genus

bash ../bac.genus.cmd.txt 2>&1 | tee log_cmd.txt

for d in `find . -mindepth 1 -maxdepth 1 -type d | sort `; do
    echo "echo \"==> Processing ${d} <==\""
    echo bash ${d}/1_pair.sh;
    echo bash ${d}/2_rawphylo.sh;
    echo bash ${d}/3_multi.sh;
    echo ;
done  > runall.sh

sh runall.sh 2>&1 | tee log_runall.txt

# Clean
```

# Summary

## Copy xlsx files

```bash
mkdir -p ~/data/bacteria/bac_summary/xlsx
cd ~/data/bacteria/bac_summary/xlsx

find  ~/data/bacteria/bac.working -type f -name "*.common.xlsx" |
    grep -v "vs[A-Z]" |
    parallel cp {} .

find  ~/data/bacteria/bac.working -type f -name "*.gc.xlsx" |
    grep -v "vs[A-Z]" |
    parallel cp {} .

```

## Genome list

Create `bac.list.csv` from `bac.WORKING.csv` with sequence lengths.

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
        s/_/ /;
        print $_;
    }
    ' \
    >> species_all.lst

echo "#target" > target_all.lst
perl -l -MPath::Tiny -e '
    BEGIN {
        @ls = map {/^#/ and s/^(#+\s*\w+).*/\1/; $_} 
            map {s/\w+,//; $_} 
            map {s/^###\s*//; $_} 
            path(q{~/Scripts/withncbi/doc/bac_target_OG.md})->lines({chomp => 1});
    }
    
    for (@ls) { 
        (/^\s*$/ or /^##\s+/ or /^#\s+(\w+)/) and next; 
        print $_;
    }
    ' \
    >> target_all.lst

echo "#abbr,species,accession,length" > length.tmp
find ~/data/bacteria/bac.working -type f -name "chr.sizes" |
    sort |
    parallel --jobs 1 --keep-order -r '
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

echo "#abbr,subgroup,genus,species,taxon_id" > abbr.tmp
cat ~/data/bacteria/bac_summary/bac.WORKING.csv |
    grep -v "^#" |
    perl -nla -F"," -e 'print qq{$F[7],$F[4],$F[3],$F[2],$F[0]}' \
    >> abbr.tmp

# #abbr,species,accession,length,subgroup,genus,species,taxon_id
cat length.tmp abbr.tmp |
    perl ~/Scripts/withncbi/util/merge_csv.pl \
        -f 0 --concat -o stdout |
    perl -nl -a -F"," -e 'print qq{$F[4],$F[5],$F[6],$F[0],$F[7],$F[2],$F[3]}' \
    > list.tmp

echo "#subgroup,genus,species,abbr,taxon_id,accession,length" > bac.list.csv
cat list.tmp |
    grep -v "#" |
    perl -nl -a -F',' -MPath::Tiny -e '
        BEGIN{
            %species, %target;
            my @l1 = path(q{species_all.lst})->lines({ chomp => 1});
            $species{$l1[$_]} = $_ for (0 .. $#l1);
            my @l2 = path(q{target_all.lst})->lines({ chomp => 1});
            $target{$l2[$_]} = $_ for (0 .. $#l2);
        }
        my $idx_s = $species{$F[2]};
        die qq{$_\n} unless defined $idx_s;
        my $idx_t = exists $target{$F[3]} ? $target{$F[3]} : 999_999;
        print qq{$_,$idx_s,$idx_t};
    ' |
    sort -n -t',' -k8,8 -k9,9 |
    cut -d',' -f 1-7 \
    >> bac.list.csv

rm *.tmp
```

## Genome alignment statistics

Some species will be filtered out here.

Criteria:

* Coverage >= 0.4
* Total number of indels >= 100
* D of multiple alignments < 0.2

```bash
mkdir -p ~/data/bacteria/bac_summary/table

cd ~/data/bacteria/bac_summary/xlsx
cat <<EOF > Table_alignment.tt
---
autofit: A:F
texts:
  - text: "Species"
    pos: A1
  - text: "No. of genomes"
    pos: B1
  - text: "Genome size on average (Mb)"
    pos: C1
  - text: "Aligned length (Mb)"
    pos: D1
  - text: "Indels per 100 bp"
    pos: E1
  - text: "Substitutions per 100 bp"
    pos: F1
  - text: "D on average"
    pos: G1
  - text: "GC-content"
    pos: H1
  - text: "Coverage on average"
    pos: I1
  - text: "Indels"
    pos: J1
[% FOREACH item IN data -%]
  - text: [% item.name %]
    pos: A[% loop.index + 2 %]
  - text: =D[% loop.index + 2 %]/C[% loop.index + 2 %]
    pos: I[% loop.index + 2 %]
  - text: =E[% loop.index + 2 %]/100*D[% loop.index + 2 %]*1000*1000
    pos: J[% loop.index + 2 %]
[% END -%]
borders:
  - range: A1:J1
    top: 1
    bottom: 1
ranges:
[% FOREACH item IN data -%]
  [% item.file %]:
    basic:
      - copy: B2
        paste: B[% loop.index + 2 %]
      - copy: B3
        paste: C[% loop.index + 2 %]
      - copy: B4
        paste: D[% loop.index + 2 %]
      - copy: B5
        paste: E[% loop.index + 2 %]
      - copy: B6
        paste: F[% loop.index + 2 %]
      - copy: B7
        paste: G[% loop.index + 2 %]
      - copy: B8
        paste: H[% loop.index + 2 %]
[% END -%]
EOF

cat ~/data/bacteria/bac_summary/table/species_all.lst \
    | grep -v "^#" \
    | TT_FILE=Table_alignment.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, { name => $_, file => qq{$species.common.xlsx}, }; 
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, })
                or die Template->error;
        }
    ' \
    > Table_alignment_all.yml

perl ~/Scripts/fig_table/xlsx_table.pl -i Table_alignment_all.yml

# Under Windows for Excel formula
perl d:/Scripts/fig_table/xlsx2xls.pl --csv -d d:/data/bacteria/bac_summary/xlsx/Table_alignment_all.xlsx

# Back to Mac
perl -pi -e 's/\r\n/\n/g;' Table_alignment_all.csv
cp -f Table_alignment_all.xlsx ~/data/bacteria/bac_summary/table
cp -f Table_alignment_all.csv ~/data/bacteria/bac_summary/table

# real filter
cd ~/data/bacteria/bac_summary/table
cat Table_alignment_all.csv \
    | perl -nla -F',' -e '
        $F[0] =~ s/"//g;
        print $F[0] if ($F[1] !~ /[\.\d]+/ or $F[8] < 0.4 or $F[9] < 100 or $F[6] > 0.2);
    ' \
    > species_exclude.lst
    
cat Table_alignment_all.csv \
    | perl -nla -F',' -e '
        $F[0] =~ s/"//g;
        print $F[0] if ($F[1] =~ /[\.\d]+/ and $F[7] >= 0.4 and $F[7] <= 0.5);
    ' \
    > species_all_demo.lst

grep -v -Fx -f species_exclude.lst species_all.lst > species.lst
grep -Fx -f species_all_demo.lst species.lst > species_demo.lst

#
cd ~/data/bacteria/bac_summary/xlsx
cat ~/data/bacteria/bac_summary/table/species.lst \
    | grep -v "^#" \
    | TT_FILE=Table_alignment.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, { name => $_, file => qq{$species.common.xlsx}, }; 
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, })
                or die Template->error;
        }
    ' \
    > Table_alignment.yml

perl ~/Scripts/fig_table/xlsx_table.pl -i Table_alignment.yml

cp -f Table_alignment.xlsx ~/data/bacteria/bac_summary/table
```

Table_S_bac for GC

```bash
cd ~/data/bacteria/bac_summary/xlsx

cat <<'EOF' > Table_S_bac.tt
---
autofit: A:J
texts:
  - text: "Species"
    pos: A2:A3
    merge: 1
  - text: "No. of genomes"
    pos: B2:B3
    merge: 1
  - text: "GC-content"
    pos: C2:C3
    merge: 1
  - text: "Aligned length (Mb)"
    pos: D2:D3
    merge: 1
  - text: "Indels per 100 bp"
    pos: E2:E3
    merge: 1
  - text: "SNPs per 100 bp"
    pos: F2:F3
    merge: 1
  - text: "Correlation coefficients (r) between"
    pos: G2:J2
    merge: 1
  - text: "GC & D"
    pos: G3
  - text: "GC & Indel"
    pos: H3
  - text: "CV & D"
    pos: I3
  - text: "CV & Indel"
    pos: J3
  - text: "r_squared"
    pos: K2:N2
    merge: 1
  - text: "GC & D"
    pos: K3
  - text: "GC & Indel"
    pos: L3
  - text: "CV & D"
    pos: M3
  - text: "CV & Indel"
    pos: N3
  - text: "p_value"
    pos: O2:R2
    merge: 1
  - text: "GC & D"
    pos: O3
  - text: "GC & Indel"
    pos: P3
  - text: "CV & D"
    pos: Q3
  - text: "CV & Indel"
    pos: R3
  - text: "slope"
    pos: S2:V2
    merge: 1
  - text: "GC & D"
    pos: S3
  - text: "GC & Indel"
    pos: T3
  - text: "CV & D"
    pos: U3
  - text: "CV & Indel"
    pos: V3
  - text: =CONCATENATE(IF(S4<0,"-",""),ROUND(SQRT(K4),3),IF(O4<0.01,"**",IF(O4<0.05,"*","")))
    pos: G4
[% FOREACH item IN data -%]
  - text: [% item.name %]
    pos: A[% loop.index + 4 %]
[% END -%]
borders:
  - range: A2:J2
    top: 1
  - range: A3:J3
    bottom: 1
ranges:
[% FOREACH item IN data -%]
  [% item.common_file %]:
    basic:
      - copy: B2
        paste: B[% loop.index + 4 %]
      - copy: B8
        paste: C[% loop.index + 4 %]
      - copy: B4
        paste: D[% loop.index + 4 %]
      - copy: B5
        paste: E[% loop.index + 4 %]
      - copy: B6
        paste: F[% loop.index + 4 %]
  [% item.gc_file %]:
    segment_gc_indel_3:
      - copy: R2
        paste: K[% loop.index + 4 %]
      - copy: R3
        paste: O[% loop.index + 4 %]
      - copy: R5
        paste: S[% loop.index + 4 %]
      - copy: R20
        paste: L[% loop.index + 4 %]
      - copy: R21
        paste: P[% loop.index + 4 %]
      - copy: R23
        paste: T[% loop.index + 4 %]
    segment_cv_indel_3:
      - copy: R2
        paste: M[% loop.index + 4 %]
      - copy: R3
        paste: Q[% loop.index + 4 %]
      - copy: R5
        paste: U[% loop.index + 4 %]
      - copy: R20
        paste: N[% loop.index + 4 %]
      - copy: R21
        paste: R[% loop.index + 4 %]
      - copy: R23
        paste: V[% loop.index + 4 %]
[% END -%]
EOF

cat ~/data/bacteria/bac_summary/table/species_all.lst \
    | grep -v "^#" \
    | TT_FILE=Table_S_bac.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, {
            name        => $_,
            common_file => qq{$species.common.xlsx},
            gc_file     => qq{$species.gc.xlsx},
        };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, })
                or die Template->error;
        }
    ' \
    > Table_S_bac_all.yml
perl ~/Scripts/fig_table/xlsx_table.pl -i Table_S_bac_all.yml

cat ~/data/bacteria/bac_summary/table/species.lst \
    | grep -v "^#" \
    | TT_FILE=Table_S_bac.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, {
            name        => $_,
            common_file => qq{$species.common.xlsx},
            gc_file     => qq{$species.gc.xlsx},
        };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, })
                or die Template->error;
        }
    ' \
    > Table_S_bac.yml
perl ~/Scripts/fig_table/xlsx_table.pl -i Table_S_bac.yml

```

NCBI Taxonomy tree

```bash
mkdir -p ~/data/bacteria/bac_summary/group
cd ~/data/bacteria/bac_summary/group

cat ~/data/bacteria/bac_summary/table/species.lst \
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
    > species_tree.sh

bash species_tree.sh > species.tree

```

## sep_chart of d1, d2

`collect_xlsx.pl`

```bash
cd ~/data/bacteria/bac_summary/xlsx

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

cat ~/data/bacteria/bac_summary/table/species.lst \
    | grep -v "^#" \
    | TT_FILE=cmd_collect_d1_d2.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, { name => $species, }; 
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, }) 
                or die Template->error;
        }
    ' \
    > cmd_collect_d1_d2.sh

cd ~/data/bacteria/bac_summary/xlsx
bash cmd_collect_d1_d2.sh

```

`sep_chart.pl`

```bash
mkdir -p ~/data/bacteria/bac_summary/fig

cd ~/data/bacteria/bac_summary/xlsx

cat <<'EOF' > cmd_chart_d1_d2.tt
perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d1.xlsx \
    -xl "Distance to indels ({italic(d)[1]})" \
    -yl "Nucleotide diversity ({italic(D)})" \
    -xr "A2:A8" -yr "B2:B8" \
    --y_min 0.0 --y_max [% y_max %] \
    -x_min 0 -x_max 5 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot -ms

perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d1_comb.xlsx \
    -xl "Distance to indels ({italic(d)[1]})" \
    -yl "Nucleotide diversity ({italic(D)})" \
    -xr "A2:A8" -yr "B2:B8" \
    --y_min 0.0 --y_max [% y_max %] \
    -x_min 0 -x_max 5 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot

perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d2.xlsx \
    -xl "Reciprocal of indel density ({italic(d)[2]})" \
    -yl "Nucleotide diversity ({italic(D)})" \
    -xr "A2:A23" -yr "B2:B23" \
    --y_min 0.0 --y_max [% y_max2 %] \
    -x_min 0 -x_max 20 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot -ms

perl ~/Scripts/fig_table/sep_chart.pl \
    -i cmd_d2_comb.xlsx \
    -xl "Reciprocal of indel density ({italic(d)[2]})" \
    -yl "Nucleotide diversity ({italic(D)})" \
    -xr "A2:A23" -yr "B2:B23" \
    --y_min 0.0 --y_max [% y_max2 %] \
    -x_min 0 -x_max 20 \
    -rb "^([% FOREACH item IN data %][% item.name %]|[% END %]NON_EXIST)$" \
    -rs "NON_EXIST" \
    --postfix [% postfix %] --style_dot

EOF

cat ~/data/bacteria/bac_summary/table/species.lst \
    | TT_FILE=cmd_chart_d1_d2.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, { name => $species, }; 
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE},
                { data => \@data,
                y_max => 0.15,
                y_max2 => 0.15,
                postfix => all, }) 
                or die Template->error;
        }
    ' \
    > cmd_chart.sh

bash cmd_chart.sh
rm ~/data/bacteria/bac_summary/xlsx/*.csv
cp ~/data/bacteria/bac_summary/xlsx/*.pdf ~/data/bacteria/bac_summary/fig

```

## CorelDRAW GC charts

```bash
cd ~/data/bacteria/bac_summary/xlsx

# Fig_S_bac_d1_gc_cv
cat <<'EOF' > Fig_S_bac_d1_gc_cv.tt
[% USE Math -%]
---
charts:
[% FOREACH item IN data -%]
  [% item.common_file %]:
    d1_comb_pi_gc_cv:
      4:
        - [% loop.index % 8 %]
        - [% Math.int(loop.index / 8) %]
[% END -%]
texts:
[% FOREACH item IN data -%]
  - text: [% item.name %]
    size: 8
    bold: 1
    italic: 1
    pos:
      - [% loop.index % 8 %]
      - [% Math.int(loop.index / 8) - 0.05 %]
[% END -%]
EOF

cat ~/data/bacteria/bac_summary/table/species.lst \
    | grep -v "^#" \
    | TT_FILE=Fig_S_bac_d1_gc_cv.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, {
            name        => $_,
            common_file => qq{$species.common.xlsx},
            gc_file     => qq{$species.gc.xlsx},
        };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, }) 
                or die Template->error;
        }
    ' \
    > Fig_S_bac_d1_gc_cv.yml

# Under Windows
perl d:/Scripts/fig_table/corel_fig.pl -i d:/data/bacteria/bac_summary/xlsx/Fig_S_bac_d1_gc_cv.yml

# Back to Mac
# Fig_S_bac_gc_indel
cat <<'EOF' > Fig_S_bac_gc_indel.tt
[% USE Math -%]
---
charts:
[% FOREACH item IN data -%]
  [% item.gc_file %]:
    segment_gc_indel_3:
      2:
        - [% loop.index % 8 %]
        - [% Math.int(loop.index / 8) %]
[% END -%]
texts:
[% FOREACH item IN data -%]
  - text: [% item.name %]
    size: 8
    bold: 1
    italic: 1
    pos:
      - [% loop.index % 8 %]
      - [% Math.int(loop.index / 8) - 0.05 %]
[% END -%]
EOF

cat ~/data/bacteria/bac_summary/table/species.lst \
    | grep -v "^#" \
    | TT_FILE=Fig_S_bac_gc_indel.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, {
            name        => $_,
            common_file => qq{$species.common.xlsx},
            gc_file     => qq{$species.gc.xlsx},
        };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, }) 
                or die Template->error;
        }
    ' \
    > Fig_S_bac_gc_indel.yml

# Under Windows
perl d:/Scripts/fig_table/corel_fig.pl -i d:/data/bacteria/bac_summary/xlsx/Fig_S_bac_gc_indel.yml

# Back to Mac
# Fig_S_bac_gc_indel_demo
cat <<'EOF' > Fig_S_bac_gc_indel_demo.tt
[% USE Math -%]
---
charts:
[% FOREACH item IN data -%]
  [% item.gc_file %]:
    segment_gc_indel_3:
      2:
        - 0
        - [% loop.index %]
    segment_cv_indel_3:
      2:
        - 1
        - [% loop.index %]
[% END -%]
texts:
[% FOREACH item IN data -%]
  - text: [% item.name %]
    size: 8
    bold: 1
    italic: 1
    pos:
      - 0
      - [% loop.index - 0.05 %]
[% END -%]
EOF

cat ~/data/bacteria/bac_summary/table/species_demo.lst \
    | grep -v "^#" \
    | sort \
    | TT_FILE=Fig_S_bac_gc_indel_demo.tt perl -MTemplate -nl -e '
        my $species = $_;
        $species =~ s/ /_/g;
        push @data, {
            name        => $_,
            common_file => qq{$species.common.xlsx},
            gc_file     => qq{$species.gc.xlsx},
        };
        END {
            $tt = Template->new;
            $tt->process($ENV{TT_FILE}, { data => \@data, }) 
                or die Template->error;
        }
    ' \
    > Fig_S_bac_gc_indel_demo.yml

# Under Windows
perl d:/Scripts/fig_table/corel_fig.pl -i d:/data/bacteria/bac_summary/xlsx/Fig_S_bac_gc_indel_demo.yml

# Back to Mac
cp -f Fig_S_bac_d1_gc_cv.cdr ~/data/bacteria/bac_summary/fig
cp -f Fig_S_bac_gc_indel.cdr ~/data/bacteria/bac_summary/fig
cp -f Fig_S_bac_gc_indel_demo.cdr ~/data/bacteria/bac_summary/fig

```

