# Processing bacterial genomes species by species

[TOC levels=1-3]: # ""

- [Processing bacterial genomes species by species](#processing-bacterial-genomes-species-by-species)
- [Init genome report database.](#init-genome-report-database)
- [Download sequences and regenerate lineage information.](#download-sequences-and-regenerate-lineage-information)
  - [Numbers for higher ranks](#numbers-for-higher-ranks)
  - [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
  - [Exclude diverged strains](#exclude-diverged-strains)
- [Prepare sequences for lastz](#prepare-sequences-for-lastz)
- [Alignments](#alignments)
  - [Create `bac_target_OG.md` for picking targets and outgroups.](#create-bac_target_ogmd-for-picking-targets-and-outgroups)
  - [Retrieve targets' assemblies](#retrieve-targets-assemblies)
  - [Create alignments plans without outgroups](#create-alignments-plans-without-outgroups)
  - [Plans for align-able targets](#plans-for-align-able-targets)
  - [Aligning w/o outgroups](#aligning-wo-outgroups)
  - [Aligning with outgroups](#aligning-with-outgroups)
  - [Self alignments](#self-alignments)
- [Summary](#summary)
  - [Copy xlsx files](#copy-xlsx-files)
  - [Genome list](#genome-list)
  - [Statistics of genome alignments](#statistics-of-genome-alignments)
  - [sep_chart of d1, d2](#sep_chart-of-d1-d2)
  - [CorelDRAW GC charts](#coreldraw-gc-charts)


# Init genome report database.

* Create databases by following steps in
  [`db/README.md`](https://github.com/wang-q/withncbi/blob/master/db/README.md#genome-reports)

* Find valid species.

  * We got **230** species.
  * 68 species have `species_code`
  * Some species have vast numbers of strains, we will exclude `NZ_*` sequences

```shell script
mkdir -p ~/data/bacteria/summary
cd ~/data/bacteria/summary

mysql -ualignDB -palignDB gr_prok -e '
    SELECT  species_id `#species_id`,
            COUNT(*) count,
            MAX(CHAR_LENGTH(code)) species_code
    FROM gr
    WHERE   1 = 1
    AND (status LIKE "%Complete%"               # complete genomes
        OR status LIKE "%Chromosome%")
    AND species NOT LIKE "%Candidatus%"         # uncertainty classification
    AND taxonomy_id != species_id               # no strain ID
    AND organism_name NOT LIKE "%,%"            # avoid commas in names
    AND (chr IS NOT NULL OR LENGTH(CHR) > 0)    # has chromosome accession
    AND (LENGTH(wgs) = 0 OR wgs IS NULL)        # avoid bad assembly
    AND species_member > 2
    AND genus IS NOT NULL
    GROUP BY species_id
    HAVING count > 2 # AND species_code > 0     # having enough and representative members
    ORDER BY subgroup, species_id
    ' |
    cut -f 1 |
    grep -v "^#" \
    > SPECIES_ID.lst
cat SPECIES_ID.lst | wc -l

```

* Expand species to strains. (Nested single quotes in bash should be '\'')

  Got **2665** strains.

```shell script
cd ~/data/bacteria/summary

cat SPECIES_ID.lst |
    parallel --keep-order -r -j 8 '
        mysql -ualignDB -palignDB gr_prok -e '\''
            SELECT  taxonomy_id `#strain_taxonomy_id`,
                    organism_name `strain`,
                    species,
                    genus,
                    subgroup,
                    `code`,
                    chr
            FROM gr
            WHERE 1 = 1
            AND (status LIKE "%Complete%"               # complete genomes
                OR status LIKE "%Chromosome%")
            AND species NOT LIKE "%Candidatus%"         # uncertainty classification
            AND taxonomy_id != species_id               # no strain ID
            AND organism_name NOT LIKE "%,%"            # avoid commas in names
            AND (chr IS NOT NULL OR LENGTH(CHR) > 0)    # has chromosome accession
            AND (LENGTH(wgs) = 0 OR wgs IS NULL)        # avoid bad assembly
            AND species_id = {}
            AND LENGTH(chr) > 6
            ORDER BY released_date                      # oldest first
        '\'' |
            grep -v "^#" |
            perl -nla -F"," -e '\''
                $F[4] =~ s/\W+/_/g;
                $F[4] =~ s/^_//g;
                $F[4] =~ s/_$//g;
                print join q{,}, @F;
            '\'' |
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
    > STRAIN.csv
cat STRAIN.csv | wc -l

```

* Create abbreviations.

* Exclude all strains of "NZ_*" in
  * Salmonella enterica
  * Escherichia coli
  * Listeria monocytogenes
  * Helicobacter pylori
  * Chlamydia trachomatis
  * Staphylococcus aureus
  * Mycobacterium tuberculosis

```shell script
cd ~/data/bacteria/summary

echo '#strain_taxonomy_id,strain,species,genus,subgroup,code,accession,abbr' > ABBR.csv.tmp
cat STRAIN.csv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "2,3,4" -s "," -m 0 --shortsub |
    sort -t',' -k5,5 -k4,4 -k3,3 -k6,6 \
    >> ABBR.csv.tmp

cat ABBR.csv.tmp |
    cut -d, -f 3 |
    uniq -c |
    sort -nr |
    head -n 10

cat ABBR.csv.tmp |
    perl -nla -F"," -e '
        if (
            $F[2] eq q{Salmonella enterica}
            or $F[2] eq q{Escherichia coli}
            or $F[2] eq q{Listeria monocytogenes}
            or $F[2] eq q{Helicobacter pylori}
            or $F[2] eq q{Chlamydia trachomatis}
            or $F[2] eq q{Staphylococcus aureus}
            or $F[2] eq q{Mycobacterium tuberculosis}
        ) {
            $F[6] =~ /^NZ_/ and next;
        }

        print;
    ' \
    > ABBR.csv

cat ABBR.csv |
    cut -d, -f 3 |
    uniq -c |
    sort -nr |
    head -n 10

```

# Download sequences and regenerate lineage information.

We don't rename sequences here, so the file has three columns. **2203** accessions.

And create `bac_ncbi.csv` with abbr names as taxon file.

```shell script
mkdir -p ~/data/bacteria/GENOMES
cd ~/data/bacteria/GENOMES

cat ../summary/ABBR.csv |
    grep -v '^#' |
    perl -nla -F"," -e 'print qq{$F[0],$F[7]}' |
    uniq |
    perl ~/Scripts/withncbi/taxon/strain_info.pl --stdin --withname --file bac_ncbi.csv

echo "#strain_name,accession,strain_taxon_id" > bac_name_acc_id.csv
cat ../summary/ABBR.csv |
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

# Some warnings about trans-splicing genes from BioPerl, just ignore them
# eutils restricts 3 connections
cat bac_name_acc_id.csv |
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
    tee bac_seq.log

# count downloaded sequences
find . -maxdepth 2 -name "*.fa" | wc -l
find . -maxdepth 1 -type d | wc -l

# failed files
find . -maxdepth 2 -type f -size -1k | grep ".fa$"

```

## Numbers for higher ranks

18 subgroups, 113 genera and 230 species.

```shell script
cd ~/data/bacteria/summary/

# count every ranks
cat ABBR.csv | sed -e '1d' | cut -d',' -f 3 | sort | uniq > species.list
cat ABBR.csv | sed -e '1d' | cut -d',' -f 4 | sort | uniq > genus.list.tmp
cat ABBR.csv | sed -e '1d' | cut -d',' -f 5 | sort | uniq > subgroup.list
wc -l subgroup.list genus.list.tmp species.list

rm *.tmp

```

## Raw phylogenetic tree by MinHash

```shell script
mkdir -p ~/data/bacteria/mash
cd ~/data/bacteria/mash

for name in $(cat ../summary/ABBR.csv | sed -e '1d' | cut -d"," -f 8 | sort); do
    2>&1 echo "==> ${name}"

    if [[ -e ${name}.msh ]]; then
        continue
    fi

    find ../GENOMES/${name} -name "*.fa" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${name}" -o ${name}
done

mkdir -p ~/data/bacteria/summary/subgroup
cd ~/data/bacteria/summary/subgroup

for subgroup in $(cat ~/data/bacteria/summary/subgroup.list); do
    echo >&2 "==> ${subgroup}"

    mash triangle -E -p 8 -l <(
            cat ../ABBR.csv |
                tsv-filter -d, --str-eq "5:${subgroup}" |
                cut -d, -f 8 |
                parallel echo "../../mash/{}.msh"
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
                pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
            tmp <- as.matrix(tmp)
            mat <- tmp[,-1]
            rownames(mat) <- tmp[,1]

            dist_mat <- as.dist(mat)
            clusters <- hclust(dist_mat, method = "ward.D2")
            tree <- as.phylo(clusters)
            write.tree(phy=tree, file="tree.nwk")

            group <- cutree(clusters, h=0.5) # k=3
            groups <- as.data.frame(group)
            groups$ids <- rownames(groups)
            rownames(groups) <- NULL
            groups <- groups[order(groups$group), ]
            write_tsv(groups, "groups.tsv")
        '

    nw_display -s -b 'visibility:hidden' -w 600 -v 30 tree.nwk |
        rsvg-convert -o ${subgroup}.png

    mv tree.nwk ${subgroup}.nwk
    mv groups.tsv ${subgroup}.groups.tsv

done

```

* Abnormal strains in Species

```shell script
cd ~/data/bacteria/summary
find ~/data/bacteria/summary/subgroup -name "*.groups.tsv" |
    sort |
    parallel -j 1 -k '
        cat {} | sed -e "1d" | xargs -I[] echo "{/.}_[]"
    ' |
    sed -e 's/.groups_/_/' |
    sed -e $'s/ /\t/' \
    > all.groups.tsv

cat species.list |
    parallel -j 1 -k '
        group=$(
            tsv-join all.groups.tsv -d 2 \
                -f <(cat ABBR.csv | grep -w {} | cut -d, -f 8) \
                -k 1 |
                cut -f 1 |
                sort |
                uniq
        )
        number=$(echo group | wc -l)
        echo -e "{}\t$number"
    '

```

## Exclude diverged strains

Check raw trees generated by the previous step.

* 391904, Bifidobacterium longum subsp. infantis ATCC 15697 = JCM 1222 = DSM 20088, 2008-11-20,
  Complete Genome,
* 1496303,Bacillus subtilis subsp. globigii,Bacillus
  subtilis,Bacillus,Firmicutes,,NZ_CP014840.1,Baci_subtilis_globigii
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
* 1243591, Salmonella enterica subsp. enterica serovar Quebec str. S-1267

* 1385755,synthetic Escherichia coli C321.deltaA,Escherichia
  coli,Escherichia,Gammaproteobacteria,,CP006698.1,Es_coli_synthetic_Escherichia_coli_C321_deltaA

```text
SELECT taxonomy_id, organism_name, released_date, status, code
FROM gr_prok.gr
WHERE species = "Gluconobacter oxydans"
and status NOT IN ('Contig', 'Scaffold')
ORDER BY released_date, status, code

```

```shell script
cd ~/data/bacteria/summary

cat ABBR.csv |
    grep -v "391904," |
    grep -v "1496303," |
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
    grep -v "372461," |
    grep -v "1243591," |
    grep -v "1385755," \
    > WORKING.csv

cat WORKING.csv |
    cut -d, -f 3 |
    uniq -c |
    sort -nr |
    head -n 10

```

# Prepare sequences for lastz

```shell script
cd ~/data/bacteria/GENOMES

find . -maxdepth 1 -mindepth 1 -type d |
    sort |
    parallel --no-run-if-empty --linebuffer -k -j 2 '
        echo >&2 "==> {}"

        if [ -e {}/chr.fasta ]; then
            echo >&2 "    {} has been processed"
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

# Alignments

## Create `bac_target_OG.md` for picking targets and outgroups.

Manually edit it then move to `~/Scripts/withncbi/doc/bac_target_OG.md`.

* Listed targets were well-curated.

* Outgroups can be changes with fewer intentions.

```shell script
cd ~/data/bacteria/summary

cat WORKING.csv |
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

# 230
cat bac_target_OG.md |
    grep -v '^#' |
    grep -E '\S+' |
    wc -l

# mv bac_target_OG.md ~/Scripts/withncbi/doc/bac_target_OG.md

```

## Retrieve targets' assemblies

**192** targets.

```shell script
cd ~/data/bacteria/summary

cat WORKING.csv |
    grep -v "^#" |
    perl -na -F"," -MPath::Tiny -e '
        BEGIN{
            %h = ();
            @ls = grep {/\S/}
                  grep {!/^#/}
                  path(q{~/Scripts/withncbi/doc/bac_target_OG.md})->lines({chomp => 1});
            for (@ls) {
                @fs = split(/,/);
                $h{$fs[1]} = 1;
            }
            undef @ls;
        }

        chomp for @F;
        if (exists $h{$F[7]}) {
            printf qq{%s\t%s\t%s\n}, $F[2], $F[7], $F[0];
        }
    ' |
    parallel --colsep '\t' --keep-order -r -j 8 '
        mysql -ualignDB -palignDB ar_refseq -e '\''
            SELECT  "{2}" `#name`,
                    ftp_path `ftp_path`,
                    species `organism`,
                    assembly_level `assembly_level`
            FROM ar
            WHERE 1 = 1
            AND taxonomy_id = {3}
            '\'' |
        grep -v "^#"
    ' |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    keep-header -- sort -k3,3 -k1,1 \
    > target.assembly.tsv

wc -l target.assembly.tsv

```

```shell script
cd ~/data/bacteria/

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f summary/target.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/target.assembly.rsync.sh

bash ASSEMBLY/target.assembly.collect.sh

tar cvfz target.rna.tar.gz \
    ASSEMBLY/target.assembly.collect.csv \
    $(find ASSEMBLY -name "*_rna_from_genomic.fna.gz")

tar cvfz target.genomic.tar.gz \
    ASSEMBLY/target.assembly.collect.csv \
    $(find ASSEMBLY -name "*_genomic.fna.gz" | grep -v "_from_" | sort)

```

## Create alignments plans without outgroups

```text
WORKING.csv
#strain_taxonomy_id,strain,species,genus,subgroup,code,accession,abbr
```

**200** species and **21** genera.

```shell script
mkdir -p ~/data/bacteria/taxon
cd ~/data/bacteria/taxon

echo -e "#Serial\tGroup\tCount\tTarget" > group_target.tsv

cat ../summary/WORKING.csv |
    grep -v "^#" |
    SERIAL=1 perl -na -F"," -MPath::Tiny -e '
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
                    printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
                    path(qq{$name})->spew(map {qq{$_\n}} @s);
                    $ENV{SERIAL}++;
                }
            }
            $name = $F[2];
            %id_of = ();
        }
        $id_of{$F[7]} = $F[0]; # same strain multiple chromosomes collapsed here

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = $h{$name};
            printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
            path(qq{$name})->spew(map {qq{$_\n}} @s);
        }' \
    >> group_target.tsv

cat ../summary/WORKING.csv |
    grep -v "^#" |
    SERIAL=501 perl -na -F"," -MPath::Tiny -e '
        BEGIN{
            our $name = q{};
            our %id_of = ();
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
                if (scalar @s > 1) {
                    printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
                    path(qq{$name})->spew(map {qq{$_\n}} @s);
                    $ENV{SERIAL}++;
                }
            }
            $name = $genus;
            %id_of = ();
        }
        $id_of{$F[7]} = $F[0]; # multiple chromosomes collapsed here

        END {
            my @s = sort {$id_of{$a} <=> $id_of{$b}} keys %id_of;
            my $t = shift @s;
            if (scalar @s > 1) {
                printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $name, scalar @s, $t;
                path(qq{$name})->spew(map {qq{$_\n}} @s);
            }
        }
    '  \
    >> group_target.tsv

```

## Plans for align-able targets

```shell script
cd ~/data/bacteria/taxon

cat ~/Scripts/withncbi/doc/bac_target_OG.md |
    grep -v "^#" |
    grep . |
    perl -nla -F"," -e 'print $F[1] ' \
    > targets.tmp

cat ../summary/all.groups.tsv |
    grep -F -w -f targets.tmp |
    perl -nla -F"\t" -e '($g, $s) = split q{_}, $F[1]; print qq{$F[0]\t${g}_${s}}' |
    tsv-summarize --group-by 1 --count |
    tsv-filter --ge 2:2 |
    cut -f 1 \
    > subgroups.tmp

cat ../summary/all.groups.tsv |
    grep -F -w -f targets.tmp |
    grep -F -w -f subgroups.tmp |
    tsv-summarize --group-by 1 --values 2 \
    > subgroups.lst.tmp

cat subgroups.lst.tmp |
    grep -v "^#" |
    SERIAL=901 perl -na -F"\t" -MPath::Tiny -e '
        chomp for @F;
        my $group = $F[0];
        my @targets = split /\|/, $F[1];

        printf qq{%s\t%s\t%s\t%s\n}, $ENV{SERIAL}, $group, scalar @targets, $targets[0];
        path(qq{$group})->spew(map {qq{$_\n}} @targets);
        $ENV{SERIAL}++;
    '  \
    >> group_target.tsv

rm *.tmp

```

## Aligning w/o outgroups

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/bacteria/ \
    wangq@202.119.37.251:data/bacteria

rsync -avP \
    ~/Scripts/withncbi/ \
    wangq@202.119.37.251:Scripts/withncbi

# rsync -avP wangq@202.119.37.251:data/bacteria/ ~/data/bacteria

```

```shell script
cd ~/data/bacteria/

# species
# HPCC restricts total jobs to 200.
#   ERROR: Number of submitted parallel jobs exceeds! Exit.
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
            --multi -o groups/species/{2} \
            --rawphylo --parallel 24 -v

#        bash groups/species/{2}/1_pair.sh
#        bash groups/species/{2}/2_rawphylo.sh
#        bash groups/species/{2}/3_multi.sh

        bsub -q mpi -n 24 -J "{2}" "
            bash groups/species/{2}/1_pair.sh
            bash groups/species/{2}/2_rawphylo.sh
            bash groups/species/{2}/3_multi.sh
        "
    '

# genus
cat taxon/group_target.tsv |
    tsv-filter -H --ge 1:501 --le 1:900 |
    sed -e '1d' | #grep "^503" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/genus/{2} \
            --rawphylo --parallel 24 -v

#        bash groups/genus/{2}/1_pair.sh
#        bash groups/genus/{2}/2_rawphylo.sh
#        bash groups/genus/{2}/3_multi.sh

        bsub -q mpi -n 24 -J "{2}" "
            bash groups/genus/{2}/1_pair.sh
            bash groups/genus/{2}/2_rawphylo.sh
            bash groups/genus/{2}/3_multi.sh
        "
    '

# subgroup
cat taxon/group_target.tsv |
    tsv-filter -H --ge 1:901 |
    sed -e '1d' | #grep "^915" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/subgroup/{2} \
            --rawphylo --parallel 24 -v

#        bash groups/subgroup/{2}/1_pair.sh
#        bash groups/subgroup/{2}/2_rawphylo.sh
#        bash groups/subgroup/{2}/3_multi.sh

        bsub -q mpi -n 24 -J "{2}" "
            bash groups/subgroup/{2}/1_pair.sh
            bash groups/subgroup/{2}/2_rawphylo.sh
            bash groups/subgroup/{2}/3_multi.sh
        "
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

# check status
echo \
    $(find groups/species -mindepth 1 -maxdepth 1 -type d | wc -l) \
    $(find groups/species -mindepth 1 -maxdepth 3 -type f -name "*.nwk.pdf" | grep -w raw -v | wc -l)

find groups/species -mindepth 1 -maxdepth 1 -type d |
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

* Review alignments and phylogenetic trees generated in `genus/` and `subgroup/`

* Add outgroups to `bac_target_OG.md` manually.

```shell script
cd ~/data/bacteria/

# species_og
cat taxon/group_target.tsv |
    tsv-filter -H --le 1:500 |
    sed -e '1d' | #grep "^163" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        outgroup=$(
            cat ~/Scripts/withncbi/doc/bac_target_OG.md |
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
            --multi -o groups/species_og/{2} \
            --outgroup ${outgroup} \
            --rawphylo --parallel 24 -v

        bsub -q mpi -n 24 -J "{2}-og" "
            bash groups/species_og/{2}/1_pair.sh
            bash groups/species_og/{2}/2_rawphylo.sh
            bash groups/species_og/{2}/3_multi.sh
            "
    '

```

## Self alignments

```shell script
cd ~/data/bacteria/

cat taxon/group_target.tsv |
    tsv-filter -H --le 1:500 |
    sed -e '1d' | # grep "^200" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --self -o groups/self/{2} \
            --circos --parallel 24 -v

        bsub -q mpi -n 24 -J "{2}-self" "
            bash groups/self/{2}/1_self.sh
            bash groups/self/{2}/3_proc.sh
            bash groups/self/{2}/4_circos.sh
            "
    '

```

# Summary

## Copy xlsx files

```shell script
mkdir -p ~/data/bacteria/summary/xlsx
cd ~/data/bacteria/summary/xlsx

find  ~/data/bacteria/bac.working -type f -name "*.common.xlsx" |
    grep -v "vs[A-Z]" |
    parallel 'cp {} .'

find  ~/data/bacteria/bac.working -type f -name "*.gc.xlsx" |
    grep -v "vs[A-Z]" |
    parallel 'cp {} .'

```

## Genome list

Create `list.csv` from `WORKING.csv` with sequence lengths.

```shell script
mkdir -p ~/data/bacteria/summary/table
cd ~/data/bacteria/summary/table

# manually adjust orders in `bac_target_OG.md`
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
cat ~/data/bacteria/summary/bac.WORKING.csv |
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

## Statistics of genome alignments

Some species will be filtered out here.

Criteria:

* Coverage >= 0.5
* Total number of indels >= 100
* D of multiple alignments < 0.2

```shell script
mkdir -p ~/data/bacteria/summary/table

cd ~/data/bacteria/summary/xlsx
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

cat ~/data/bacteria/summary/table/species_all.lst \
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
perl d:/Scripts/fig_table/xlsx2xls.pl --csv -d d:/data/bacteria/summary/xlsx/Table_alignment_all.xlsx

# Back to Mac
perl -pi -e 's/\r\n/\n/g;' Table_alignment_all.csv
cp -f Table_alignment_all.xlsx ~/data/bacteria/summary/table
cp -f Table_alignment_all.csv ~/data/bacteria/summary/table

# real filter
cd ~/data/bacteria/summary/table
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
cd ~/data/bacteria/summary/xlsx
cat ~/data/bacteria/summary/table/species.lst \
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

cp -f Table_alignment.xlsx ~/data/bacteria/summary/table
```

Table_S_bac for GC

```shell script
cd ~/data/bacteria/summary/xlsx

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

cat ~/data/bacteria/summary/table/species_all.lst \
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

cat ~/data/bacteria/summary/table/species.lst \
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

```shell script
mkdir -p ~/data/bacteria/summary/group
cd ~/data/bacteria/summary/group

cat ~/data/bacteria/summary/table/species.lst \
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

```shell script
cd ~/data/bacteria/summary/xlsx

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

cat ~/data/bacteria/summary/table/species.lst \
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

cd ~/data/bacteria/summary/xlsx
bash cmd_collect_d1_d2.sh

```

`sep_chart.pl`

```shell script
mkdir -p ~/data/bacteria/summary/fig

cd ~/data/bacteria/summary/xlsx

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

cat ~/data/bacteria/summary/table/species.lst \
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
rm ~/data/bacteria/summary/xlsx/*.csv
cp ~/data/bacteria/summary/xlsx/*.pdf ~/data/bacteria/summary/fig

```

## CorelDRAW GC charts

```shell script
cd ~/data/bacteria/summary/xlsx

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

cat ~/data/bacteria/summary/table/species.lst \
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
perl d:/Scripts/fig_table/corel_fig.pl -i d:/data/bacteria/summary/xlsx/Fig_S_bac_d1_gc_cv.yml

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

cat ~/data/bacteria/summary/table/species.lst \
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
perl d:/Scripts/fig_table/corel_fig.pl -i d:/data/bacteria/summary/xlsx/Fig_S_bac_gc_indel.yml

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

cat ~/data/bacteria/summary/table/species_demo.lst \
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
perl d:/Scripts/fig_table/corel_fig.pl -i d:/data/bacteria/summary/xlsx/Fig_S_bac_gc_indel_demo.yml

# Back to Mac
cp -f Fig_S_bac_d1_gc_cv.cdr ~/data/bacteria/summary/fig
cp -f Fig_S_bac_gc_indel.cdr ~/data/bacteria/summary/fig
cp -f Fig_S_bac_gc_indel_demo.cdr ~/data/bacteria/summary/fig

```

