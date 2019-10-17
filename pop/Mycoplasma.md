# *Mycoplasma*


[TOC levels=1-3]: # " "
- [*Mycoplasma*](#mycoplasma)
- [Strain info](#strain-info)
- [Mycoplasma: assembly](#mycoplasma-assembly)
- [Count strains](#count-strains)
- [Mycoplasma: prepare](#mycoplasma-prepare)
- [Mycoplasma: run](#mycoplasma-run)


# Strain info

| Group                     | Species          | Species ID | Comments    | Strains |
|:--------------------------|:-----------------|-----------:|:------------|--------:|
| Mycoplasma                |                  |            |             |         |
|                           | M. genitalium    |       2097 | 生殖支原体    |       5 |
|                           | M. gallisepticum |       2096 | 鸡毒支原体    |      13 |
|                           | M. pneumoniae    |       2104 | 肺炎支原体    |      17 |
|                           | M. bovis         |      28903 | 牛支原体      |       6 |
|                           | M. canis         |      29555 | 犬支原体      |       5 |
|                           | M. fermentans    |       2115 | 发酵支原体    |       5 |
|                           | M. hyopneumoniae |       2099 | 猪肺炎支原体  |       5 |
|                           | M. hyorhinis     |       2100 | 猪鼻支原体    |       6 |
|                           | M. ovipneumoniae |      29562 | 绵羊肺炎支原体 |       4 |
| Ureaplasma                |                  |            |             |         |
|                           | M. parvum        |     134821 | 微小脲原体    |       6 |
|                           | M. urealyticum   |       2130 | 解脲脲原体    |      14 |
| Mycoplasma mycoides group |                  |            |             |         |
|                           | M. capricolum    |       2095 | 山羊支原体    |       7 |
|                           | M. leachii       |       2105 |             |       3 |
|                           | M. mycoides      |       2102 | 丝状支原体    |      13 |

* [Mycoplasma](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2093)
* [Ureaplasma](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=2129)

# Mycoplasma: assembly

```bash
export RANK_NAME=Mycoplasma

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND genus_id in (2093, 2129)
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND genus_id in (2093, 2129)
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
    > ${RANK_NAME}.assembly.tsv

# comment out unneeded assembly levels

# find potential duplicated strains or assemblies
cat ${RANK_NAME}.assembly.tsv |
    cut -f 1 |
    sort |
    uniq -c |
    sort -nr

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim ${RANK_NAME}.assembly.tsv
# cp ${RANK_NAME}.assembly.tsv ~/Scripts/withncbi/pop

# Cleaning
rm raw*.*sv

unset RANK_NAME

```

```bash
cd ~/data/alignment/Mycoplasma

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Mycoplasma.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/Mycoplasma.assembly.rsync.sh

bash ASSEMBLY/Mycoplasma.assembly.collect.sh

```

# Count strains

```bash
cd ~/data/alignment/Mycoplasma

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/Mycoplasma.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep {} |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)
    
    n_strains=$(cat ASSEMBLY/Mycoplasma.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep {} |
        cut -d" " -f 1,2 |
        sort |
        wc -l)
    
    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: "Mycoplasma" \
          "Mycoplasma genitalium" "Mycoplasma gallisepticum" "Mycoplasma pneumoniae" \
          "Mycoplasma bovis" "Mycoplasma canis" "Mycoplasma fermentans" \
          "Mycoplasma hyopneumoniae" "Mycoplasma hyorhinis" "Mycoplasma ovipneumoniae" \
          "Ureaplasma" \
          "Ureaplasma parvum" "Ureaplasma urealyticum" \
          "Mycoplasma capricolum" "Mycoplasma leachii" "Mycoplasma mycoides"

#Mycoplasma      68      157
#Mycoplasma genitalium   1       5
#Mycoplasma gallisepticum        1       13
#Mycoplasma pneumoniae   1       17
#Mycoplasma bovis        1       6
#Mycoplasma canis        1       5
#Mycoplasma fermentans   1       5
#Mycoplasma hyopneumoniae        1       6
#Mycoplasma hyorhinis    1       6
#Mycoplasma ovipneumoniae        1       4
#Ureaplasma      4       22
#Ureaplasma parvum       1       6
#Ureaplasma urealyticum  1       14
#Mycoplasma capricolum   1       7
#Mycoplasma leachii      1       3
#Mycoplasma mycoides     1       12

mkdir -p taxon

parallel --no-run-if-empty --linebuffer -k -j 4 '
    cat ASSEMBLY/Mycoplasma.assembly.collect.csv |
        cut -d"," -f 1,2 |
        grep {} |
        cut -d"," -f 1 \
        > taxon/{= $_ =~ s/ /_/g =} # replace spaces with underscore
    ' ::: "Mycoplasma genitalium" "Mycoplasma gallisepticum" "Mycoplasma pneumoniae" \
          "Mycoplasma bovis" "Mycoplasma canis" "Mycoplasma fermentans" \
          "Mycoplasma hyopneumoniae" "Mycoplasma hyorhinis" "Mycoplasma ovipneumoniae" \
          "Ureaplasma parvum" "Ureaplasma urealyticum" \
          "Mycoplasma capricolum" "Mycoplasma leachii" "Mycoplasma mycoides"

cat ASSEMBLY/Mycoplasma.assembly.collect.csv |
    tsv-filter -H -d"," --istr-ne 12:"Contig" |
    tsv-filter -H -d"," --istr-ne 12:"Scaffold" |
    tsv-filter -H -d"," --not-blank 18 |
    sed '1d' | # remove header
    cut -d"," -f 1 \
    > taxon/Mycoplasma

cat \
    taxon/Mycoplasma_genitalium \
    taxon/Mycoplasma_pneumoniae \
    > taxon/gen_pne

cat \
    taxon/Ureaplasma_parvum \
    taxon/Ureaplasma_urealyticum \
    > taxon/Ureaplasma

cat \
    taxon/Mycoplasma_capricolum \
    taxon/Mycoplasma_leachii \
    taxon/Mycoplasma_mycoides \
    > taxon/mycoides_group

wc -l taxon/*
#  22 taxon/gen_pne
#  22 taxon/mycoides_group
#  30 taxon/Mycoplasma
#   6 taxon/Mycoplasma_bovis
#   5 taxon/Mycoplasma_canis
#   7 taxon/Mycoplasma_capricolum
#   5 taxon/Mycoplasma_fermentans
#  13 taxon/Mycoplasma_gallisepticum
#   5 taxon/Mycoplasma_genitalium
#   6 taxon/Mycoplasma_hyopneumoniae
#   6 taxon/Mycoplasma_hyorhinis
#   3 taxon/Mycoplasma_leachii
#  12 taxon/Mycoplasma_mycoides
#   4 taxon/Mycoplasma_ovipneumoniae
#  17 taxon/Mycoplasma_pneumoniae
#  20 taxon/Ureaplasma
#   6 taxon/Ureaplasma_parvum
#  14 taxon/Ureaplasma_urealyticum

find taxon -maxdepth 1 -type f -not -name "*.replace.tsv" |
    sort |
    xargs -i basename {} \
    > species.list

```

# Mycoplasma: prepare

* `--perseq` for RefSeq_category Reference Genome assemblies.

```bash
cd ~/data/alignment/Mycoplasma

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    --perseq M_bovis_PG45 \
    --perseq M_canis_PG_14 \
    --perseq M_cap_capricolum_ATCC_27343 \
    --perseq M_fer_M64 \
    --perseq M_gallis_R_low \
    --perseq M_gen_G37 \
    --perseq M_hyop_J \
    --perseq M_hyor_SK76 \
    --perseq M_lea_PG50 \
    --perseq M_myc_mycoides_SC_PG1 \
    --perseq M_ovip_NM2010 \
    --perseq M_pne_M129 \
    --perseq U_par_3_ATCC_27815 \
    --perseq U_ure_10_ATCC_33699 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--parallel 12"

bash GENOMES/0_prep.sh

# gff
for n in M_bovis_PG45 M_canis_PG_14 M_cap_capricolum_ATCC_27343 \
         M_fer_M64 M_gallis_R_low M_gen_G37 \
         M_hyop_J M_hyor_SK76 M_lea_PG50 \
         M_myc_mycoides_SC_PG1 M_ovip_NM2010 M_pne_M129 \
         U_par_3_ATCC_27815 U_ure_10_ATCC_33699 ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}/")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"
    
    gzip -dcf ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

# Mycoplasma: run

```bash
# species and targets
ARRAY=(
    'gen_pne::M_gen_G37'
    'mycoides_group::M_myc_mycoides_SC_PG1'
    'Mycoplasma::M_gen_G37'
    'Mycoplasma_bovis::M_bovis_PG45'
    'Mycoplasma_canis::M_canis_PG_14'
    'Mycoplasma_capricolum::M_cap_capricolum_ATCC_27343'
    'Mycoplasma_fermentans::M_fer_M64'
    'Mycoplasma_gallisepticum::M_gallis_R_low'
    'Mycoplasma_genitalium::M_gen_G37'
    'Mycoplasma_hyopneumoniae::M_hyop_J'
    'Mycoplasma_hyorhinis::M_hyor_SK76'
    'Mycoplasma_leachii::M_lea_PG50'
    'Mycoplasma_mycoides::M_myc_mycoides_SC_PG1'
    'Mycoplasma_ovipneumoniae::M_ovip_NM2010'
    'Mycoplasma_pneumoniae::M_pne_M129'
    'Ureaplasma::U_ure_10_ATCC_33699'
    'Ureaplasma_parvum::U_par_3_ATCC_27815'
    'Ureaplasma_urealyticum::U_ure_10_ATCC_33699'
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

