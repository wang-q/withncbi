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

| #tax_id | species                     | RS  | comments                |
|---------|-----------------------------|-----|-------------------------|
|         |                             |     | M. tuberculosis complex |
| 78331   | M. canettii                 | 29  | 卡氏分枝杆菌                  |
| 1773    | M. tuberculosis             | 286 | 结核分枝杆菌                  |
|         |                             |     | M. avium complex        |
| 1764    | M. avium                    | 39  | 鸟分枝杆菌                   |
| 339268  | M. colombiense              | 9   |                         |
| 1767    | M. intracellulare           | 62  | 胞内分枝杆菌                  |
| 560555  | M. mantenii                 | 3   |                         |
| 701042  | M. marseillense             | 5   |                         |
| 1138383 | M. paraintracellulare       | 12  | 副胞内分枝杆菌                 |
| 701043  | M. timonense                | 1   |                         |
|         |                             |     | M. simiae complex       |
| 1964395 | M. ahvazicum                | 1   |                         |
| 761804  | M. europaeum                | 2   |                         |
| 292462  | M. florentinum              | 2   |                         |
| 36812   | M. genavense                | 1   |                         |
| 53376   | M. heidelbergense           | 2   |                         |
| 33895   | M. interjectum              | 1   |                         |
| 120959  | M. kubicae                  | 6   |                         |
| 141349  | M. lentiflavum              | 1   |                         |
| 154654  | M. montefiorense            | 1   |                         |
| 767916  | M. paraense                 | 4   |                         |
| 240125  | M. parascrofulaceum         | 1   |                         |
| 185642  | M. parmense                 | 2   |                         |
| 220927  | M. saskatchewanense         | 2   |                         |
| 243061  | M. sherrisii                | 1   |                         |
| 722731  | M. shigaense                | 3   |                         |
| 1784    | M. simiae                   | 3   | 猿分支杆菌                   |
| 470076  | M. stomatepiae              | 1   |                         |
| 47839   | M. triplex                  | 2   |                         |
|         |                             |     | Others                  |
| 1927124 | M. aquaticum                | 1   |                         |
| 1790    | M. asiaticum                | 4   | 亚洲分枝杆菌                  |
| 2094119 | M. basiliense               | 1   |                         |
| 56425   | M. bohemicum                | 2   | 波希米亚分枝杆菌                |
| 84962   | M. botniense                | 1   |                         |
| 1273442 | M. bourgelatii              | 1   |                         |
| 43348   | M. branderi                 | 2   |                         |
| 28045   | M. celatum                  | 2   |                         |
| 44010   | M. conspicuum               | 2   |                         |
| 1775    | M. cookii                   | 1   |                         |
| 1430326 | M. decipiens                | 1   |                         |
| 1801    | M. diernhoferi              | 3   |                         |
| 482462  | M. dioxanotrophicus         | 1   |                         |
| 1260918 | M. fragae                   | 2   |                         |
| 117567  | M. frederiksbergense        | 1   |                         |
| 1778    | M. gordonae                 | 1   | 戈登分枝杆菌                  |
| 1552759 | M. grossiae                 | 1   |                         |
| 29311   | M. haemophilum              | 5   |                         |
| 110505  | M. heckeshornense           | 3   |                         |
| 512402  | M. heraklionense            | 1   |                         |
| 49897   | M. hodleri                  | 2   |                         |
| 152142  | M. holsaticum               | 1   |                         |
| 1768    | M. kansasii                 | 17  | 堪萨斯分枝杆菌                 |
| 1069220 | M. koreense                 | 2   |                         |
| 2212479 | M. kyogaense                | 1   |                         |
| 169765  | M. lacus                    | 1   |                         |
| 2048550 | M. lehmannii                | 2   |                         |
| 1769    | M. leprae                   | 3   | 麻风分枝杆菌                  |
| 480418  | M. lepromatosis             | 1   |                         |
| 261524  | M. liflandii                | 1   |                         |
| 1780    | M. malmoense                | 2   | 莫尔门分枝杆菌                 |
| 1781    | M. marinum                  | 19  | 海洋分枝杆菌                  |
| 244292  | M. nebraskense              | 1   |                         |
| 242737  | M. neglectum                | 1   |                         |
| 2048551 | M. neumannii                | 1   |                         |
| 459858  | M. noviomagense             | 1   |                         |
| 2492438 | M. novum                    | 1   |                         |
| 1841861 | M. numidiamassiliense       | 1   |                         |
| 2738409 | M. ostraviense              | 2   |                         |
| 370524  | M. pallens                  | 1   |                         |
| 53378   | M. paraffinicum             | 1   |                         |
| 1389713 | M. paragordonae             | 6   |                         |
| 590652  | M. paraseoulense            | 2   |                         |
| 1487726 | M. persicum                 | 4   |                         |
| 2341080 | M. pseudokansasii           | 1   |                         |
| 265949  | M. pseudoshottsii           | 1   |                         |
| 1841860 | M. rhizamassiliense         | 1   |                         |
| 486698  | M. riyadhense               | 9   |                         |
| 1796    | M. senegalense              | 4   |                         |
| 386911  | M. seoulense                | 1   |                         |
| 29313   | M. shimoidei                | 2   |                         |
| 398694  | M. shinjukuense             | 1   |                         |
| 133549  | M. shottsii                 | 2   |                         |
| 627089  | M. simulans                 | 2   |                         |
| 1785    | M. sp.                      | 58  |                         |
| 886343  | M. spongiae                 | 1   |                         |
| 1908205 | M. syngnathidarum           | 2   |                         |
| 1841859 | M. terramassiliense         | 1   |                         |
| 2162698 | M. uberis                   | 1   |                         |
| 1809    | M. ulcerans                 | 3   |                         |
| 1719132 | M. vicinigordonae           | 1   |                         |
| 547163  | M. vulneris                 | 1   |                         |
| 1789    | M. xenopi                   | 3   | 蟾分枝杆菌                   |
|         |                             |     | Actinobacteria          |
| 33910   | Amycolatopsis mediterranei  | 4   | 地中海拟无枝酸菌                |
| 1717    | Corynebacterium diphtheriae | 16  | 白喉棒状杆菌                  |
| 1912    | Streptomyces hygroscopicus  | 3   | 吸水链霉菌                   |

| #tax_id | biotype                            | RS  |            |
|---------|------------------------------------|-----|:-----------|
| 33894   | M. tuberculosis variant africanum  | 33  | 非洲分枝杆菌     |
| 1765    | M. tuberculosis variant bovis      | 102 | 牛分枝杆菌, 卡介苗 |
| 115862  | M. tuberculosis variant caprae     | 3   | 山羊分枝杆菌     |
| 1806    | M. tuberculosis variant microti    | 5   | 仓鼠分枝杆菌     |
| 194542  | M. tuberculosis variant pinnipedii | 1   | 鳍足分枝杆菌     |

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
    tsv-filter -H --ge RS:3 --or --ge GB:3 |
    sed 's/Mycobacterium /M. /g' |
    mlr --itsv --omd cat

cat species.count.tsv |
    tsv-filter -H --lt RS:1 --ge GB:2 |
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

| #tax_id | species     | RS  | GB  |
|---------|-------------|-----|-----|
| 134601  | M. goodii   | 0   | 7   |
| 39693   | M. porcinum | 0   | 6   |
| 318424  | M. rufum    | 0   | 3   |

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

| #tax_id | group                               | RS   |
|---------|-------------------------------------|------|
| 1844474 | M. mungi                            | 1    |
| 1305738 | M. orygis                           | 2    |
| 78331   | M. canettii                         | 34   |
| 1773    | M. tuberculosis                     | 6830 |
| 194542  | M. tuberculosis variant pinnipedii  | 3    |
| 182785  | M. tuberculosis subsp. tuberculosis | 1    |
| 115862  | M. tuberculosis variant caprae      | 4    |
| 33894   | M. tuberculosis variant africanum   | 34   |
| 1806    | M. tuberculosis variant microti     | 6    |
| 1765    | M. tuberculosis variant bovis       | 143  |

M. tuberculosis var. microti is distinct from other M. tuberculosis biotypes

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
        AND genus IN ('Amycolatopsis', 'Corynebacterium', 'Streptomyces')
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

| #tax_id | species                            | RS  |
|---------|------------------------------------|-----|
| 33910   | Amycolatopsis mediterranei         | 4   |
| 1717    | Corynebacterium diphtheriae        | 14  |
| 1718    | Corynebacterium glutamicum         | 9   |
| 1719    | Corynebacterium pseudotuberculosis | 14  |
| 65058   | Corynebacterium ulcerans           | 5   |
| 43771   | Corynebacterium urealyticum        | 2   |
| 29303   | Streptomyces cattleya              | 2   |
| 1969    | Streptomyces chartreusis           | 2   |
| 1912    | Streptomyces hygroscopicus         | 2   |
| 1916    | Streptomyces lividans              | 2   |
| 1927    | Streptomyces rimosus               | 2   |
| 54571   | Streptomyces venezuelae            | 2   |

## MTBC: assembly

```shell
cd ~/data/alignment/MTBC

# Not M. tuberculosis, M. avium and M. intracellulare
# Good M. sp.
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND genus IN ('Mycobacterium')
        AND (
            (species LIKE '% sp.%' AND assembly_level IN ('Complete Genome', 'Chromosome') )
            OR
            (species NOT LIKE '% sp.%')
        )
        AND species != 'Mycobacterium tuberculosis'
        AND species != 'Mycobacterium avium'
        AND species != 'Mycobacterium intracellulare'
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > raw.tsv

# M. avium and M. intracellulare
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species IN ('Mycobacterium avium', 'Mycobacterium intracellulare')
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
        AND (
            (organism_name LIKE '% subsp %')
            OR
            (organism_name NOT LIKE '% variant %' AND assembly_level IN ('Complete Genome', 'Chromosome') )
            OR
            (organism_name LIKE '% variant %'
                AND (
                    (organism_name NOT LIKE '% variant bovis %')
                    OR
                    (assembly_level IN ('Complete Genome', 'Chromosome') )
                )
            )
        )
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    >> raw.tsv

# Missing in refseq
#134601 | M. goodii
#39693 | M. porcinum
#318424 | M. rufum
echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN (134601, 39693, 318424)
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite \
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

# Comment out
# M_col_Mycobacterium_tuberculosis_TKK_01_0051_GCF_000661085_1

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

#find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | head |
#    tr "/" "\t" |
#    cut -f 2 |
#    tsv-join --exclude -k 1 -f ASSEMBLY/rsync.tsv -d 1 |
#    xargs -I[] rm -fr ASSEMBLY/[]

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
#  710 ASSEMBLY/n50.pass.csv
#  965 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/MTBC.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/MTBC.assembly.pass.csv

wc -l ASSEMBLY/MTBC.assembly*csv
#   965 ASSEMBLY/MTBC.assembly.collect.csv
#   710 ASSEMBLY/MTBC.assembly.pass.csv

```

```shell
cd ~/data/alignment/MTBC

# Group by genus
cat ASSEMBLY/MTBC.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 2 |
    tsv-select -d" " -f 1 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            tsv-select -d, -f 2 |
            grep {} |
            tsv-select -d" " -f 1,2 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            tsv-select -d, -f 2 |
            grep {} |
            tsv-select -d" " -f 1,2 |
            sort |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    '
#Amycolatopsis   1       4
#Corynebacterium 1       16
#Mycobacterium   90      686
#Streptomyces    1       3

# Group by species
cat ASSEMBLY/MTBC.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 2 |
    tsv-select -d" " -f 1,2 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_strains=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            tsv-select -d, -f 2 |
            grep {} |
            tsv-select -d" " -f 1,2 |
            sort |
            wc -l)

        printf "%s\t%d\n" {} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 4,1,2 |
    nwr append stdin -r "species group" |
    sed 's/NA//g' |
    tsv-sort -k4,4 -k2,2 |
    sed 's/Mycobacterium /M. /g' |
    (echo -e '#tax_id\tspecies\tRS\t' && cat) |
    mlr --itsv --omd cat
# Paste this table to the top of this document

cat ASSEMBLY/MTBC.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 2 |
    grep "Mycobacterium tuberculosis variant" |
    cut -d" " -f 1-4 |
    tsv-uniq |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_strains=$(cat ASSEMBLY/MTBC.assembly.pass.csv |
            tsv-select -d, -f 2 |
            grep {} |
            tsv-select -d" " -f 1,2 |
            sort |
            wc -l)

        printf "%s\t%d\n" {} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 4,1,2 |
    sed 's/Mycobacterium /M. /g' |
    (echo -e '#tax_id\tbiotype\tRS' && cat) |
    mlr --itsv --omd cat
# Paste this table to the top of this document

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

