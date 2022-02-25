# *Mycobacterium tuberculosis* complex

[TOC levels=1-3]: # ""

- [*Mycobacterium tuberculosis* complex](#mycobacterium-tuberculosis-complex)
  - [Strain info](#strain-info)
  - [MTBC: assembly](#mtbc-assembly)
  - [Count strains](#count-strains)
  - [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
  - [NCBI taxonomy](#ncbi-taxonomy)
  - [Groups and targets](#groups-and-targets)
  - [MTBC: prepare](#mtbc-prepare)
  - [MTBC: run](#mtbc-run)

## Strain info

* [Mycobacterium](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1763)
* [Mycobacterium tuberculosis complex](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=77643)

| Group                   | Species               |  Tax ID | Comments |
|:------------------------|:----------------------|--------:|:---------|
| M. tuberculosis complex |                       |         |          |
|                         | M. canetti            |   78331 | 卡氏分枝杆菌   |
|                         | M. decipiens          | 1430326 |          |
|                         | M. mungi              | 1844474 |          |
|                         | M. orygis             | 1305738 |          |
|                         | M. tuberculosis       |    1773 | 结核分枝杆菌   |
| M. avium complex        |                       |         |          |
|                         | M. avium              |    1764 | 鸟分枝杆菌    |
|                         | M. chimaera           |  222805 |          |
|                         | M. colombiense        |  339268 |          |
|                         | M. intracellulare     |    1767 | 胞内分枝杆菌   |
|                         | M. mantenii           |  560555 |          |
|                         | M. marseillense       |  701042 |          |
|                         | M. paraintracellulare | 1138383 | 副胞内分枝杆菌  |
| M. simiae complex       |                       |         |          |
|                         | M. genavense          |   36812 |          |
|                         | M. parascrofulaceum   |  240125 |          |
|                         | M. simiae             |    1784 | 猿分支杆菌    |
| Close to MSC            |                       |         |          |
|                         | M. asiaticum          |    1790 | 亚洲分枝杆菌   |
|                         | M. bohemicum          |    1998 | 波希米亚分枝杆菌 |
| Others                  |                       |         |          |
|                         | M. gordonae           |    1778 | 戈登分枝杆菌   |
|                         | M. haemophilum        |   29311 |          |
|                         | M. kansasii           |    1955 | 堪萨斯分枝杆菌  |
|                         | M. leprae             |    1769 | 麻风分枝杆菌   |
|                         | M. liflandii          |  261524 |          |
|                         | M. malmoense          |    1780 | 莫尔门分枝杆菌  |
|                         | M. marinum            |    1781 | 海洋分枝杆菌   |
|                         | M. persicum           | 1487726 |          |
|                         | M. pseudoshottsii     |  265949 |          |
|                         | M. ulcerans           |    1809 |          |
|                         | M. xenopi             |    1789 | 蟾分枝杆菌    |

### List all ranks

```shell
mkdir -p ~/data/alignment/MTBC
cd ~/data/alignment/MTBC

nwr member Mycobacterium |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Mycobacterium -r "species group" -r subspecies -r biotype |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3r |
    sed 's/Mycobacterium /M. /g' |
    mlr --itsv --omd cat

```

| rank          | count |
|---------------|-------|
| genus         | 1     |
| species       | 149   |
| no rank       | 39    |
| species group | 3     |
| strain        | 3002  |
| subspecies    | 9     |
| biotype       | 5     |

| #tax_id | sci_name                                | rank          |
|---------|-----------------------------------------|---------------|
| 1124626 | M. ulcerans subsp. shinshuense          | subspecies    |
| 1203599 | M. intracellulare subsp. yongonense     | subspecies    |
| 1770    | M. avium subsp. paratuberculosis        | subspecies    |
| 182785  | M. tuberculosis subsp. tuberculosis     | subspecies    |
| 222805  | M. intracellulare subsp. chimaera       | subspecies    |
| 35617   | M. intracellulare subsp. intracellulare | subspecies    |
| 439334  | M. avium subsp. hominissuis             | subspecies    |
| 44282   | M. avium subsp. silvaticum              | subspecies    |
| 44454   | M. avium subsp. avium                   | subspecies    |
| 120793  | M. avium complex (MAC)                  | species group |
| 2249310 | M. simiae complex                       | species group |
| 77643   | M. tuberculosis complex                 | species group |
| 115862  | M. tuberculosis variant caprae          | biotype       |
| 1765    | M. tuberculosis variant bovis           | biotype       |
| 1806    | M. tuberculosis variant microti         | biotype       |
| 194542  | M. tuberculosis variant pinnipedii      | biotype       |
| 33894   | M. tuberculosis variant africanum       | biotype       |

### Species with assemblies

```shell
cd ~/data/alignment/MTBC

SPECIES=$(
    nwr member Mycobacterium -r species |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        sed '1d' |
        cut -f 1 |
        sort
)

for S in $SPECIES; do
    RS=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    GB=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
            " |
            sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
    )

    if [[ ${RS} -gt 0 ]] || [[ ${GB} -gt 0 ]]; then
        echo -e "$S\t$RS\t$GB"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr |
    (echo -e '#tax_id\tspecies\tRS\tGB' && cat) \
    > species.count.tsv

cat species.count.tsv |
    tsv-filter -H --ge RS:5 --or --ge GB:5 |
    sed 's/Mycobacterium /M. /g' |
    mlr --itsv --omd cat

```

| #tax_id | species               | RS   | GB   |
|---------|-----------------------|------|------|
| 1773    | M. tuberculosis       | 6830 | 7074 |
| 1764    | M. avium              | 247  | 417  |
| 1767    | M. intracellulare     | 83   | 86   |
| 1768    | M. kansasii           | 34   | 47   |
| 78331   | M. canettii           | 34   | 35   |
| 1781    | M. marinum            | 26   | 26   |
| 339268  | M. colombiense        | 20   | 20   |
| 1138383 | M. paraintracellulare | 14   | 15   |
| 1487726 | M. persicum           | 10   | 12   |
| 1790    | M. asiaticum          | 10   | 10   |
| 486698  | M. riyadhense         | 10   | 10   |
| 1780    | M. malmoense          | 8    | 8    |
| 1809    | M. ulcerans           | 7    | 11   |
| 1778    | M. gordonae           | 7    | 9    |
| 120959  | M. kubicae            | 7    | 7    |
| 1389713 | M. paragordonae       | 6    | 6    |
| 560555  | M. mantenii           | 6    | 6    |
| 701042  | M. marseillense       | 6    | 6    |
| 1769    | M. leprae             | 5    | 6    |
| 29311   | M. haemophilum        | 5    | 5    |
| 1789    | M. xenopi             | 4    | 6    |
| 512402  | M. heraklionense      | 1    | 7    |
| 547163  | M. vulneris           | 1    | 5    |
| 134601  | M. goodii             | 0    | 7    |
| 39693   | M. porcinum           | 0    | 6    |

### MTBC assemblies

```shell
cd ~/data/alignment/MTBC

GROUP=$(
    nwr member "Mycobacterium tuberculosis complex" -r species -r subspecies -r biotype |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        sed '1d' |
        cut -f 1
)

for G in $GROUP; do
    printf "$G\t"
    echo "
        SELECT
            tax_id
        FROM ar
        WHERE 1=1
            AND genus_id = 1763
        " |
        sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
        parallel -j 4 "
            IS_ANCESTOR=\$(nwr lineage {} | tsv-filter --str-eq 3:$G)
            if [[ \${#IS_ANCESTOR} -gt 1 ]]; then
                echo {} | nwr append stdin
            fi
            " |
        wc -l
done |
    nwr append stdin |
    tsv-select -f 1,3,2 |
    sed 's/Mycobacterium /M. /g' |
    (echo -e '#tax_id\tgroup\tRS' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | group                               | RS   | comments   |
|---------|-------------------------------------|------|:-----------|
| 1844474 | M. mungi                            | 1    |            |
| 1305738 | M. orygis                           | 2    |            |
| 78331   | M. canettii                         | 34   | 卡氏分枝杆菌     |
| 1773    | M. tuberculosis                     | 6830 |            |
| 194542  | M. tuberculosis variant pinnipedii  | 3    | 鳍足分枝杆菌     |
| 182785  | M. tuberculosis subsp. tuberculosis | 1    |            |
| 115862  | M. tuberculosis variant caprae      | 4    | 山羊分枝杆菌     |
| 33894   | M. tuberculosis variant africanum   | 34   | 非洲分枝杆菌     |
| 1806    | M. tuberculosis variant microti     | 6    | 仓鼠分枝杆菌     |
| 1765    | M. tuberculosis variant bovis       | 143  | 牛分枝杆菌, 卡介苗 |

M. tuberculosis var. microti is distinct from other M. tuberculosis strains

### MTB biotypes

56 strains and 190 assemblies have biotype information.

```shell
cd ~/data/alignment/MTBC

nwr member "Mycobacterium tuberculosis"

nwr member "Mycobacterium tuberculosis" -r strain |
    sed '1d' | # head -n 1000 |
    cut -f 1 |
    parallel -j 4 '
        IS_ANCESTOR=$(nwr lineage {} | tsv-filter --str-eq 1:biotype)
        if [[ ${#IS_ANCESTOR} -gt 1 ]]; then
            echo {} | nwr append stdin
        fi
        ' |
    wc -l
# 56

echo "
    SELECT
        tax_id
    FROM ar
    WHERE 1=1
        AND species_id = 1773
        AND tax_id != 1773
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    parallel -j 4 '
        IS_ANCESTOR=$(nwr lineage {} | tsv-filter --str-eq 1:biotype)
        if [[ ${#IS_ANCESTOR} -gt 1 ]]; then
            echo {} | nwr append stdin
        fi
        ' |
    wc -l
# 190

```

### Other species of Actinobacteria

```shell
cd ~/data/alignment/MTBC

echo "
    SELECT
        tax_id
    FROM ar
    WHERE 1=1
        AND genus IN ('Amycolatopsis', 'Corynebacterium', 'Streptomyces') -- species_id IN (33910, 1717, 1912)
        AND assembly_level IN ('Complete Genome', 'Chromosome')
        AND tax_id != species_id
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    parallel -j 4 '
        IS_STRAIN=$(nwr info --tsv {} | tsv-filter --str-eq 3:strain)
        if [[ ${#IS_STRAIN} -gt 1 ]]; then
            echo {} | nwr append stdin -r species --id
        fi
        ' |
    tsv-summarize -g 3,2 --count |
    tsv-filter --gt 3:1 |
    tsv-sort -k2,2 |
    (echo -e '#tax_id\tspecies\tRS' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | species                            | RS  | comments |
|---------|------------------------------------|-----|:---------|
| 33910   | Amycolatopsis mediterranei         | 4   | 地中海拟无枝酸菌 |
| 1717    | Corynebacterium diphtheriae        | 14  | 白喉棒状杆菌   |
| 1718    | Corynebacterium glutamicum         | 9   |          |
| 1719    | Corynebacterium pseudotuberculosis | 14  |          |
| 65058   | Corynebacterium ulcerans           | 5   |          |
| 43771   | Corynebacterium urealyticum        | 2   |          |
| 29303   | Streptomyces cattleya              | 2   |          |
| 1969    | Streptomyces chartreusis           | 2   |          |
| 1912    | Streptomyces hygroscopicus         | 2   | 吸水链霉菌    |
| 1916    | Streptomyces lividans              | 2   |          |
| 1927    | Streptomyces rimosus               | 2   |          |
| 54571   | Streptomyces venezuelae            | 2   |          |

## MTBC: assembly

```shell
cd ~/data/alignment/MTBC

# Not M. tuberculosis and M. avium
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus IN ('Mycobacterium')
        AND species != 'Mycobacterium tuberculosis'
        AND species != 'Mycobacterium avium'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > raw.tsv

# M. avium
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species = 'Mycobacterium avium'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# M. tuberculosis biotypes
#1773    | M. tuberculosis
#194542  | M. tuberculosis variant pinnipedii
#182785  | M. tuberculosis subsp. tuberculosis
#115862  | M. tuberculosis variant caprae
#33894   | M. tuberculosis variant africanum
#1806    | M. tuberculosis variant microti
#1765    | M. tuberculosis variant bovis
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species = 'Mycobacterium tuberculosis'
        AND tax_id != 1773
        AND ( (
            organism_name LIKE '% variant %'
            ) OR (
            organism_name LIKE '% subsp %'
            ) OR (
            organism_name NOT LIKE '% variant %' AND assembly_level IN ('Complete Genome', 'Chromosome')
            ) )
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# Good outgroups
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN (33910, 1717, 1912)
        AND assembly_level IN ('Complete Genome', 'Chromosome')
        AND tax_id != species_id
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

cat raw.tsv |
    grep -v '^#' |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[5]}++;
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > MTBC.assembly.tsv

# find potential duplicated strains or assemblies
cat MTBC.assembly.tsv |
    tsv-uniq -f 1 --repeated

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim MTBC.assembly.tsv
# cp MTBC.assembly.tsv ~/Scripts/withncbi/pop

# Cleaning
rm raw*.*sv

```

```shell
cd ~/data/alignment/MTBC

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/MTBC.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/MTBC.assembly.rsync.sh

bash ASSEMBLY/MTBC.assembly.collect.sh

```

## Count strains

```shell
cd ~/data/alignment/MTBC

for dir in $(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | sort); do
    1>&2 echo "==> ${dir}"
    name=$(basename ${dir})

    find ${dir} -type f -name "*_genomic.fna.gz" |
        grep -v "_from_" | # exclude CDS and rna
        xargs cat |
        faops n50 -C -S stdin |
        (echo -e "name\t${name}" && cat) |
        datamash transpose
done |
    tsv-uniq |
    tee ASSEMBLY/n50.tsv

cat ASSEMBLY/n50.tsv |
    tsv-filter \
        -H --or \
        --le 4:100 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50*
#     264 ASSEMBLY/n50.pass.csv
#     286 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/MTBC.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/MTBC.assembly.pass.csv

wc -l ASSEMBLY/MTBC.assembly*csv
#     286 ASSEMBLY/MTBC.assembly.collect.csv
#     264 ASSEMBLY/MTBC.assembly.pass.csv

```

```shell
cd ~/data/alignment/MTBC

cat ASSEMBLY/MTBC.assembly.pass.csv |
    sed -e '1d' |
    cut -d"," -f 2 |
    grep -v "Candidatus" |
    cut -d" " -f 1 |
    sort |
    uniq |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            grep "{}" |
            cut -d" " -f 1,2 |
            sort |
            uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            grep "{}" |
            cut -d" " -f 1,2 |
            sort |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    '

#Amycolatopsis	1	3
#Corynebacterium	1	24
#Mycobacterium	18	232
#Streptomyces	1	4

cat ASSEMBLY/MTBC.assembly.pass.csv |
    sed -e '1d' |
    cut -d"," -f 2 |
    grep -v "Candidatus" |
    cut -d" " -f 1,2 |
    sort |
    uniq |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            grep {} |
            cut -d" " -f 1,2 |
            sort |
            uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            grep {} |
            cut -d" " -f 1,2 |
            sort |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    '
#Amycolatopsis mediterranei      1       3
#Corynebacterium diphtheriae     1       24
#Mycobacterium asiaticum 1       1
#Mycobacterium avium     1       7
#Mycobacterium bohemicum 1       1
#Mycobacterium canettii  1       5
#Mycobacterium colombiense       1       1
#Mycobacterium genavense 1       1
#Mycobacterium haemophilum       1       1
#Mycobacterium intracellulare    1       8
#Mycobacterium kansasii  1       4
#Mycobacterium leprae    1       3
#Mycobacterium liflandii 1       1
#Mycobacterium marinum   1       4
#Mycobacterium parascrofulaceum  1       1
#Mycobacterium pseudoshottsii    1       1
#Mycobacterium simiae    1       1
#Mycobacterium tuberculosis      1       187
#Mycobacterium ulcerans  1       2
#Mycobacterium xenopi    1       3
#Streptomyces hygroscopicus      1       4

cat ASSEMBLY/MTBC.assembly.pass.csv |
    sed -e '1d' |
    cut -d"," -f 2 |
    grep -v "Candidatus" |
    grep "Mycobacterium tuberculosis variant" |
    cut -d" " -f 1-4 |
    sort |
    uniq |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            grep {} |
            cut -d" " -f 1,2 |
            sort |
            uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            grep {} |
            cut -d" " -f 1,2 |
            sort |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    '

#Mycobacterium tuberculosis variant africanum    1       26
#Mycobacterium tuberculosis variant bovis        1       14
#Mycobacterium tuberculosis variant caprae       1       1
#Mycobacterium tuberculosis variant microti      1       1
#Mycobacterium tuberculosis variant pinnipedii   1       1

```

## Raw phylogenetic tree by MinHash

```shell
mkdir -p ~/data/alignment/MTBC/mash
cd ~/data/alignment/MTBC/mash

for name in $(cat ../ASSEMBLY/MTBC.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 ); do
    2>&1 echo "==> ${name}"

    if [[ -e ${name}.msh ]]; then
        continue
    fi

    find ../ASSEMBLY/${name} -name "*.fsa_nt.gz" -or -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${name}" -o ${name}
done

mash triangle -E -p 8 -l <(
    cat ../ASSEMBLY/MTBC.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 | parallel echo "{}.msh"
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

        group <- cutree(clusters, h=0.4) # k=5
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        write_tsv(groups, "groups.tsv")
    '

nw_display -s -b 'visibility:hidden' -w 600 -v 30 tree.nwk |
    rsvg-convert -o ~/Scripts/withncbi/image/MTBC.png

```

## NCBI taxonomy

```shell
cd ~/data/alignment/MTBC/taxon

bp_taxonomy2tree.pl -e \
    $(
        cat ASSEMBLY/MTBC.assembly.pass.csv |
            sed -e '1d' |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            cut -d" " -f 1,2 |
            tr " " "_" |
            sort |
            uniq |
            parallel echo '-s {}'
    ) \
    -s Mycobacterium_canettii \
    -s Mycobacterium_decipiens \
    -s Mycobacterium_mungi \
    -s Mycobacterium_orygis \
    -s Mycobacterium_tuberculosis \
    -s Mycobacterium_avium \
    -s Mycobacterium_chimaera \
    -s Mycobacterium_colombiense \
    -s Mycobacterium_intracellulare \
    -s Mycobacterium_mantenii \
    -s Mycobacterium_marseillense \
    -s Mycobacterium_paraintracellulare \
    -s Mycobacterium_asiaticum \
    -s Mycobacterium_gordonae \
    -s Mycobacterium_haemophilum \
    -s Mycobacterium_kansasii \
    -s Mycobacterium_leprae \
    -s Mycobacterium_malmoense \
    -s Mycobacterium_marseillense \
    -s Mycobacterium_persicum |
    sed 's/Mycobacterium/M/g' |
    sed 's/ (MAC)//g' \
    > ncbi.nwk

nw_display -s -b 'visibility:hidden' -w 600 -v 30 ncbi.nwk |
    rsvg-convert -o ~/Scripts/withncbi/image/Mycobacterium.ncbi.png

```

![Mycobacterium.ncbi.png](../image/Mycobacterium.ncbi.png)

![Mycobacterium.raxml.png](../image/Mycobacterium.raxml.png)

## Groups and targets

Review `ASSEMBLY/MTBC.assembly.pass.csv` and `mash/groups.tsv`

| #Serial | Group                                        | Count | Target                           |
|:--------|:---------------------------------------------|------:|:---------------------------------|
| 1       | A_med                                        |     3 | A_med_U32                        |
| 2       | C_dip                                        |    24 | C_dip_NCTC_13129                 |
| 3       | M_avi_int                                    |    22 | M_avi_paratuberculosis_K_10      |
| 4       | M_can_tub                                    |   191 | M_tub_H37Rv                      |
| 5       | M_hae_lep                                    |     4 | M_lep_TN                         |
| 6       | M_kan_xen                                    |     7 | M_kan_ATCC_12478                 |
| 7       | M_lif_mar_pse_ulc                            |     8 | M_mar_M                          |
| 8       | S_hyg                                        |     4 | S_hyg_jinggangensis_5008         |
| 101     | Mycobacterium                                |    20 | M_tub_H37Rv                      |
| 102     | Mycobacterium_avium                          |     7 | M_avi_paratuberculosis_K_10      |
| 103     | Mycobacterium_canettii                       |     5 | M_can_CIPT_140010059             |
| 104     | Mycobacterium_intracellulare                 |     8 | M_int_ATCC_13950                 |
| 105     | Mycobacterium_kansasii                       |     4 | M_kan_ATCC_12478                 |
| 106     | Mycobacterium_leprae                         |     3 | M_lep_TN                         |
| 107     | Mycobacterium_marinum                        |     4 | M_mar_E11                        |
| 108     | Mycobacterium_tuberculosis                   |   156 | M_tub_H37Rv                      |
| 109     | Mycobacterium_tuberculosis_variant_africanum |    26 | M_tub_variant_africanum_GM041182 |
| 110     | Mycobacterium_tuberculosis_variant_bovis     |    14 | M_tub_variant_bovis_AF2122_97    |

```shell
mkdir -p ~/data/alignment/MTBC/taxon
cd ~/data/alignment/MTBC/taxon

cp ../mash/tree.nwk .
cp ../mash/groups.tsv .

## manually removes some assemblies
#cat ../mash/groups.tsv |
#    grep -v "Ba_mic" |
#    grep -v "Cry_bai" |
#    grep -v "En_inv_IP1" |
#    grep -v "L_sp_A" \
#    > groups.tsv
#echo -e "2\tBa_mic" >> groups.tsv
#echo -e "2\tBa_mic_RI" >> groups.tsv

ARRAY=(
    'A_med::A_med_U32'
    'C_dip::C_dip_NCTC_13129'
    'M_avi_int::M_avi_paratuberculosis_K_10'
    'M_can_tub::M_tub_H37Rv'
    'M_hae_lep::M_lep_TN'
    'M_kan_xen::M_kan_ATCC_12478'
    'M_lif_mar_pse_ulc::M_mar_M'
    'S_hyg::S_hyg_jinggangensis_5008'
)

echo -e "#Serial\tGroup\tCount\tTarget" > group_target.tsv

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$(
        cat ../mash/groups.tsv |
            tsv-filter --str-eq 2:${TARGET_NAME} |
            tsv-select -f 1
    )

    cat groups.tsv |
        tsv-filter --str-eq 1:${SERIAL} |
        tsv-select -f 2 \
        > ${GROUP_NAME}

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> group_target.tsv

done

# Custom groups
ARRAY=(
    'Mycobacterium::M_tub_H37Rv'
    'Mycobacterium_avium::M_avi_paratuberculosis_K_10'
    'Mycobacterium_canettii::M_can_CIPT_140010059'
    'Mycobacterium_intracellulare::M_int_ATCC_13950'
    'Mycobacterium_kansasii::M_kan_ATCC_12478'
    'Mycobacterium_leprae::M_lep_TN'
    'Mycobacterium_marinum::M_mar_E11'
    'Mycobacterium_tuberculosis::M_tub_H37Rv'
    'Mycobacterium_tuberculosis_variant_africanum::M_tub_variant_africanum_GM041182'
    'Mycobacterium_tuberculosis_variant_bovis::M_tub_variant_bovis_AF2122_97'
)

SERIAL=100
for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$((SERIAL + 1))
    GROUP_NAME2=$(echo $GROUP_NAME | tr "_" " ")

    if [ "$GROUP_NAME" = "Mycobacterium" ]; then
        cat ../ASSEMBLY/MTBC.assembly.pass.csv |
            tsv-filter -H -d"," --not-blank 18 |
            sed '1d' |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
        echo "A_med_U32" >> ${GROUP_NAME}
        echo "C_dip_NCTC_13129" >> ${GROUP_NAME}
        echo "S_hyg_jinggangensis_5008" >> ${GROUP_NAME}
        echo "M_tub_variant_microti_OV254" >> ${GROUP_NAME}
    elif [ "$GROUP_NAME" = "Mycobacterium_tuberculosis" ]; then
        cat ../ASSEMBLY/MTBC.assembly.pass.csv |
            tsv-filter -H -d"," --istr-ne 12:"Contig" |
            tsv-filter -H -d"," --istr-ne 12:"Scaffold" |
            cut -d"," -f 1,2 |
            grep "${GROUP_NAME2}" |
            grep -v "M_tub_variant_microti" |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
    else
        cat ../ASSEMBLY/MTBC.assembly.pass.csv |
            cut -d"," -f 1,2 |
            grep "${GROUP_NAME2}" |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
    fi

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> group_target.tsv

done

mlr --itsv --omd cat group_target.tsv

cat <<'EOF' > chr-level.list
M_avi_paratuberculosis_K_10
M_can_CIPT_140010059
M_col_CECT_3035
M_hae_DSM_44634
M_int_ATCC_13950
M_int_intracellulare_MTCC_9506
M_int_yongonense_05_1390
M_kan_ATCC_12478
M_mar_E11
M_tub_H37Rv
M_tub_variant_africanum_GM041182
M_tub_variant_bovis_AF2122_97
EOF

```

## MTBC: prepare

* Rsync to hpcc

```shell
rsync -avP \
    ~/data/alignment/MTBC/ \
    wangq@202.119.37.251:data/alignment/MTBC

# rsync -avP wangq@202.119.37.251:data/alignment/MTBC/ ~/data/alignment/MTBC

```

* `--perseq` for Chromosome-level assemblies and targets

```shell
cd ~/data/alignment/MTBC/

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    $( cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 | parallel -j 1 echo " --perseq {} " ) \
    $( cat taxon/chr-level.list | parallel -j 1 echo " --perseq {} " ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--parallel 24"

bsub -q mpi -n 24 -J "MTBC-0_prep" "bash GENOMES/0_prep.sh"

ls -t output.* | head -n 1 | xargs tail -f | grep "==>"

# gff
for n in $(cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 ) \
    $( cat taxon/chr-level.list ) \
    ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"

    gzip -d -c ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

## MTBC: run

```shell
cd ~/data/alignment/MTBC/

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/{2}/ \
            --tree taxon/tree.nwk \
            --parallel 24 -v

        bsub -q mpi -n 24 -J "{2}-1_pair" "bash groups/{2}/1_pair.sh"
        bsub -w "ended({2}-1_pair)" \
            -q mpi -n 24 -J "{2}-3_multi" "bash groups/{2}/3_multi.sh"
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

```

* C_dip_NCTC_13129 as outgroup

```shell
cd ~/data/alignment/MTBC/

egaz template \
    GENOMES/M_tub_H37Rv \
    $(cat taxon/Mycobacterium | grep -v -x "M_tub_H37Rv" | xargs -I[] echo "GENOMES/[]") \
    --multi -o groups/Mycobacterium/ \
    --multiname OG --outgroup C_dip_NCTC_13129 \
    --tree groups/Mycobacterium/Results/Mycobacterium.nwk \
    --parallel 16 -v

bash groups/Mycobacterium/3_multi.sh

nw_display -s -b 'visibility:hidden' -w 600 -v 30 groups/Mycobacterium/Results/OG.nwk |
    rsvg-convert -o ~/Scripts/withncbi/image/Mycobacterium.raxml.png

```

