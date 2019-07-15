# *Mycobacterium tuberculosis* complex


[TOC levels=1-3]: # " "
- [*Mycobacterium tuberculosis* complex](#mycobacterium-tuberculosis-complex)
- [Strain info](#strain-info)
- [MTBC: assembly](#mtbc-assembly)
- [NCBI taxonomy](#ncbi-taxonomy)
- [Count strains](#count-strains)
- [MTBC: prepare](#mtbc-prepare)
- [MTBC: run](#mtbc-run)


# Strain info

| Group                   | Species               | Species ID | Comments    | Strains |
|:------------------------|:----------------------|-----------:|:------------|--------:|
| M. tuberculosis complex |                       |            |             |         |
|                         | M. canetti            |      78331 | 卡氏分枝杆菌  |       9 |
|                         | M. decipiens          |    1430326 |             |         |
|                         | M. mungi              |    1844474 |             |         |
|                         | M. orygis             |    1305738 |             |       1 |
|                         | M. tuberculosis       |       1773 | 结核分枝杆菌  |     196 |
| M. avium complex        |                       |            |             |         |
|                         | M. avium              |       1764 | 鸟分枝杆菌    |      56 |
|                         | M. chimaera           |     222805 |             |         |
|                         | M. colombiense        |     339268 |             |       1 |
|                         | M. intracellulare     |       1767 | 胞内分枝杆菌  |       9 |
|                         | M. mantenii           |     560555 |             |         |
|                         | M. marseillense       |     701042 |             |         |
|                         | M. paraintracellulare |    1138383 |             |         |
| Others                  |                       |            |             |         |
|                         | M. asiaticum          |       1790 | 亚洲分枝杆菌  |       1 |
|                         | M. gordonae           |       1778 | 戈登分枝杆菌  |         |
|                         | M. haemophilum        |      29311 |             |       1 |
|                         | M. kansasii           |       1955 | 堪萨斯分枝杆菌 |       4 |
|                         | M. leprae             |       1769 | 麻风分枝杆菌  |       5 |
|                         | M. malmoense          |       1780 | 莫尔门分枝杆菌 |         |
|                         | M. marinum            |       1781 | 海洋分枝杆菌  |       4 |
|                         | M. marseillense       |     701042 |             |         |
|                         | M. persicum           |    1487726 |             |         |

| Subspecies                      |     ID | Comments        | Strains |
|:--------------------------------|-------:|:----------------|--------:|
| M. tuberculosis var. africanum  |  33894 | 非洲分枝杆菌      |      26 |
| M. tuberculosis var. bovis      |   1765 | 牛分枝杆菌, 卡介苗 |      22 |
| M. tuberculosis var. caprae     | 115862 |                 |       1 |
| M. tuberculosis var. microti    |   1806 | 仓鼠分枝杆菌      |       2 |
| M. tuberculosis var. pinnipedii | 194542 |                 |       1 |

* [Mycobacterium](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1763)
* [Mycobacterium tuberculosis complex](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=77643)

# MTBC: assembly

```bash
export RANK_NAME=MTBC

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND ( genus_id in (1763) OR species_id in (78331, 1844474, 1305738, 1773) )
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND ( genus_id in (1763) OR species_id in (78331, 1844474, 1305738, 1773) )
    " \
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
    > raw.0.assembly.tsv

# find potential duplicated strains or assemblies
cat raw.0.assembly.tsv |
    cut -f 1 |
    sort |
    uniq -c |
    sort -nr

cat raw.0.assembly.tsv |
    tsv-filter -H --str-ne 3:"Mycobacterium tuberculosis" \
    > raw.1.assembly.tsv

cat raw.0.assembly.tsv |
    tsv-filter -H --str-eq 3:"Mycobacterium tuberculosis" |
    tsv-filter -H --str-in-fld 1:"_variant_" \
    > raw.2.assembly.tsv

cat raw.0.assembly.tsv |
    tsv-filter -H --str-eq 3:"Mycobacterium tuberculosis" |
    tsv-filter -H --istr-ne 4:"Contig" |
    tsv-filter -H --istr-ne 4:"Scaffold" \
    > raw.3.assembly.tsv

tsv-append -H raw.{1,2,3}.assembly.tsv |
    tsv-uniq -H > ${RANK_NAME}.assembly.tsv

# comment out unneeded assembly levels

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim ${RANK_NAME}.assembly.tsv
# cp ${RANK_NAME}.assembly.tsv ~/Scripts/withncbi/pop

# Cleaning
rm raw*.*sv

unset RANK_NAME

```

```bash
cd ~/data/alignment/MTBC

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/MTBC.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/MTBC.assembly.rsync.sh

bash ASSEMBLY/MTBC.assembly.collect.sh

```

# NCBI taxonomy

```bash
cd ~/data/alignment/MTBC

bp_taxonomy2tree.pl -e \
    -s "Mycobacterium canettii" \
    -s "Mycobacterium decipiens" \
    -s "Mycobacterium mungi" \
    -s "Mycobacterium orygis" \
    -s "Mycobacterium tuberculosis" \
    -s "Mycobacterium avium" \
    -s "Mycobacterium chimaera" \
    -s "Mycobacterium colombiense" \
    -s "Mycobacterium intracellulare" \
    -s "Mycobacterium mantenii" \
    -s "Mycobacterium marseillense" \
    -s "Mycobacterium paraintracellulare" \
    -s "Mycobacterium asiaticum" \
    -s "Mycobacterium gordonae" \
    -s "Mycobacterium haemophilum" \
    -s "Mycobacterium kansasii" \
    -s "Mycobacterium leprae" \
    -s "Mycobacterium malmoense" \
    -s "Mycobacterium malmoense" \
    -s "Mycobacterium marseillense" \
    -s "Mycobacterium persicum" |
    sed 's/ (MAC)//' \
    > Mycobacterium.newick

nw_display -w 600 -s Mycobacterium.newick |
    rsvg-convert -o ~/Scripts/withncbi/image/Mycobacterium.png

```

![Mycobacterium.png](../image/Mycobacterium.png)

# Count strains

```bash
cd ~/data/alignment/MTBC

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/MTBC.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)
    
    n_strains=$(cat ASSEMBLY/MTBC.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        wc -l)
    
    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: "Mycobacterium" \
          "Mycobacterium canettii" "Mycobacterium orygis" "Mycobacterium tuberculosis" \
          "Mycobacterium tuberculosis variant africanum" \
          "Mycobacterium tuberculosis variant bovis" \
          "Mycobacterium tuberculosis variant caprae" \
          "Mycobacterium tuberculosis variant microti" \
          "Mycobacterium tuberculosis variant pinnipedii" \
          "Mycobacterium avium" "Mycobacterium colombiense" \
          "Mycobacterium intracellulare" \
          "Mycobacterium asiaticum" "Mycobacterium haemophilum" \
          "Mycobacterium kansasii" "Mycobacterium leprae" \
          "Mycobacterium marinum"

#Mycobacterium   23      303
#Mycobacterium canettii  1       9
#Mycobacterium orygis    1       1
#Mycobacterium tuberculosis      1       196
#Mycobacterium tuberculosis variant africanum    1       26
#Mycobacterium tuberculosis variant bovis        1       22
#Mycobacterium tuberculosis variant caprae       1       1
#Mycobacterium tuberculosis variant microti      1       2
#Mycobacterium tuberculosis variant pinnipedii   1       1
#Mycobacterium avium     1       56
#Mycobacterium colombiense       1       1
#Mycobacterium intracellulare    1       9
#Mycobacterium asiaticum 1       1
#Mycobacterium haemophilum       1       1
#Mycobacterium kansasii  1       4
#Mycobacterium leprae    1       5
#Mycobacterium marinum   1       4

mkdir -p taxon

parallel --no-run-if-empty --linebuffer -k -j 4 '
    cat ASSEMBLY/MTBC.assembly.collect.csv |
        cut -d"," -f 1,2 |
        grep "{}" |
        cut -d"," -f 1 \
        > taxon/{= $_ =~ s/ /_/g =} # replace spaces with underscore
    ' ::: "Mycobacterium canettii" \
          "Mycobacterium tuberculosis variant africanum" \
          "Mycobacterium tuberculosis variant bovis" \
          "Mycobacterium avium" \
          "Mycobacterium intracellulare" \
          "Mycobacterium kansasii" \
          "Mycobacterium leprae" \
          "Mycobacterium marinum"

parallel --no-run-if-empty --linebuffer -k -j 4 '
    cat ASSEMBLY/MTBC.assembly.collect.csv |
        tsv-filter -H -d"," --istr-ne 12:"Contig" |
        tsv-filter -H -d"," --istr-ne 12:"Scaffold" |
        cut -d"," -f 1,2 |
        grep "{}" |
        cut -d"," -f 1 \
        > taxon/{= $_ =~ s/ /_/g =} # replace spaces with underscore
    ' ::: "Mycobacterium tuberculosis"

cat ASSEMBLY/MTBC.assembly.collect.csv |
    tsv-filter -H -d"," --not-blank 18 |
    sed '1d' | # remove header
    cut -d"," -f 1 \
    > taxon/Mycobacterium

wc -l taxon/*
#  20 taxon/Mycobacterium
#  56 taxon/Mycobacterium_avium
#   9 taxon/Mycobacterium_canettii
#   9 taxon/Mycobacterium_intracellulare
#   4 taxon/Mycobacterium_kansasii
#   5 taxon/Mycobacterium_leprae
#   4 taxon/Mycobacterium_marinum
# 155 taxon/Mycobacterium_tuberculosis
#  26 taxon/Mycobacterium_tuberculosis_variant_africanum
#  22 taxon/Mycobacterium_tuberculosis_variant_bovis

find taxon -maxdepth 1 -type f -not -name "*.replace.tsv" |
    sort |
    xargs -i basename {} \
    > species.list

```

# MTBC: prepare

* `--perseq` for RefSeq_category Reference Genome assemblies.

```bash
cd ~/data/alignment/MTBC

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    --perseq M_avi_paratuberculosis_K_10 \
    --perseq M_can_CIPT_140010059 \
    --perseq M_int_ATCC_13950 \
    --perseq M_kan_ATCC_12478 \
    --perseq M_lep_TN \
    --perseq M_mar_E11 \
    --perseq M_tub_H37Rv \
    --perseq M_tub_variant_africanum_GM041182 \
    --perseq M_tub_variant_bovis_AF2122_97 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--parallel 12"

bash GENOMES/0_prep.sh

# gff
for n in M_avi_paratuberculosis_K_10 M_can_CIPT_140010059 \
         M_int_ATCC_13950 M_kan_ATCC_12478 \
         M_lep_TN M_mar_E11 \
         M_tub_H37Rv \
         M_tub_variant_africanum_GM041182 M_tub_variant_bovis_AF2122_97; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}/")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"
    
    gzip -dcf ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

# MTBC: run

```bash
# species and targets
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

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    echo "==> ${GROUP_NAME}"
    
    egaz template \
        GENOMES/${TARGET_NAME} \
        $(cat taxon/${GROUP_NAME} | grep -v "${TARGET_NAME}" | parallel -r echo "GENOMES/{}") \
        --multi -o species/${GROUP_NAME}/ \
        --rawphylo --parallel 12 -v

    bash species/${GROUP_NAME}/1_pair.sh
    bash species/${GROUP_NAME}/2_rawphylo.sh
    bash species/${GROUP_NAME}/3_multi.sh

done

# clean
find species -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find species -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

