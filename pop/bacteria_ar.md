# Aligning various genera missing from `bacteria_gr.md`


[TOC levels=1-3]: # " "
- [Aligning various genera missing from `bacteria_gr.md`](#aligning-various-genera-missing-from-bacteria_grmd)
- [Strain info](#strain-info)
- [bacteria_ar: assembly](#bacteria_ar-assembly)
- [Count strains](#count-strains)
- [bacteria_ar: prepare](#bacteria_ar-prepare)
- [bacteria_ar: run](#bacteria_ar-run)


# Strain info

| Group               | Species                     | Species ID | Comments      | Strains |
|:--------------------|:----------------------------|-----------:|:--------------|--------:|
| Gammaproteobacteria |                             |            |               |         |
|                     | Acinetobacter junii         |      40215 | 琼氏不动杆菌    |       6 |
|                     | Citrobacter freundii        |        546 | 弗劳地枸橼酸杆菌 |       9 |
|                     | Klebsiella aerogenes        |        548 | 产气肠杆菌      |      15 |
|                     | Klebsiella oxytoca          |        571 | 产酸克雷伯菌    |      16 |
|                     | Morganella morganii         |        582 | 摩根摩根菌      |       6 |
|                     | Proteus mirabilis           |        584 | 奇异变形杆菌    |       7 |
|                     | Serratia marcescens         |        615 | 粘质沙雷菌      |      26 |
| Firmicutes/Bacilli  |                             |            |               |         |
|                     | Staphylococcus capitis      |      29388 | 头状葡萄球菌    |       7 |
|                     | Staphylococcus haemolyticus |       1283 | 溶血葡萄球菌    |       3 |
|                     | Staphylococcus hominis      |       1290 | 人葡萄球菌      |       6 |

* Synonym: Klebsiella aerogenes ==> Enterobacter aerogenes

# bacteria_ar: assembly

```shell script
export RANK_NAME=bacteria_ar

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND species_id in (40215, 546, 548, 571, 582, 584, 615, 29388, 1283, 1290)
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND species_id in (40215, 546, 548, 571, 582, 584, 615, 29388, 1283, 1290)
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
cd ~/data/alignment/bacteria_ar

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/bacteria_ar.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/bacteria_ar.assembly.rsync.sh

bash ASSEMBLY/bacteria_ar.assembly.collect.sh

```

# Count strains

```bash
cd ~/data/alignment/bacteria_ar

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/bacteria_ar.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)

    n_strains=$(cat ASSEMBLY/bacteria_ar.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        wc -l)

    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: "Acinetobacter junii" "Citrobacter freundii" "Klebsiella aerogenes" "Klebsiella oxytoca" \
          "Morganella morganii" "Proteus mirabilis" "Serratia marcescens" \
          "Staphylococcus capitis" "Staphylococcus haemolyticus" "Staphylococcus hominis"

#Acinetobacter junii     1       6
#Citrobacter freundii    1       9
#Klebsiella aerogenes    1       15
#Klebsiella oxytoca      1       16
#Morganella morganii     1       6
#Proteus mirabilis       1       7
#Serratia marcescens     1       26
#Staphylococcus capitis  1       7
#Staphylococcus haemolyticus     1       3
#Staphylococcus hominis  1       6

mkdir -p taxon

parallel --no-run-if-empty --linebuffer -k -j 4 '
    cat ASSEMBLY/bacteria_ar.assembly.collect.csv |
        cut -d"," -f 1,2 |
        grep "{}" |
        cut -d"," -f 1 \
        > taxon/{= $_ =~ s/ /_/g =} # replace spaces with underscore
    ' ::: "Acinetobacter junii" "Citrobacter freundii" "Klebsiella aerogenes" "Klebsiella oxytoca" \
          "Morganella morganii" "Proteus mirabilis" "Serratia marcescens" \
          "Staphylococcus capitis" "Staphylococcus haemolyticus" "Staphylococcus hominis"

wc -l taxon/*

find taxon -maxdepth 1 -type f -not -name "*.replace.tsv" |
    sort |
    xargs -i basename {} \
    > species.list

```

# bacteria_ar: prepare

* `--perseq` for RefSeq_category Reference Genome assemblies.

```bash
cd ~/data/alignment/bacteria_ar

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    --perseq A_jun_SH205 \
    --perseq C_fre_CFNIH1 \
    --perseq K_aer_KCTC_2190 \
    --perseq K_oxy_KONIH1 \
    --perseq M_mor_morganii_KT \
    --perseq P_mir_HI4320 \
    --perseq Se_mar_marcescens_Db11 \
    --perseq St_cap_capitis \
    --perseq St_hae_JCSC1435 \
    --perseq St_hom_hominis_C80 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--parallel 12"

bash GENOMES/0_prep.sh

# gff
for n in A_jun_SH205 C_fre_CFNIH1 K_aer_KCTC_2190 K_oxy_KONIH1 \
         M_mor_morganii_KT P_mir_HI4320 Se_mar_marcescens_Db11 \
         St_cap_capitis St_hae_JCSC1435 St_hom_hominis_C80; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"

    gzip -dcf ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

# bacteria_ar: run

```bash
# species and targets
ARRAY=(
    'Acinetobacter_junii::A_jun_SH205'
    'Citrobacter_freundii::C_fre_CFNIH1'
    'Klebsiella_aerogenes::K_aer_KCTC_2190'
    'Klebsiella_oxytoca::K_oxy_KONIH1'
    'Morganella_morganii::M_mor_morganii_KT'
    'Proteus_mirabilis::P_mir_HI4320'
    'Serratia_marcescens::Se_mar_marcescens_Db11'
    'Staphylococcus_capitis::St_cap_capitis'
    'Staphylococcus_haemolyticus::St_hae_JCSC1435'
    'Staphylococcus_hominis::St_hom_hominis_C80'
)

for item in "${ARRAY[@]}" ; do
    SPECIES_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    echo "==> ${SPECIES_NAME}"

    egaz template \
        GENOMES/${TARGET_NAME} \
        $(cat taxon/${SPECIES_NAME} | grep -v "${TARGET_NAME}" | parallel -r echo "GENOMES/{}") \
        --multi -o species/${SPECIES_NAME}/ \
        --rawphylo --parallel 12 -v

    bash species/${SPECIES_NAME}/1_pair.sh
    bash species/${SPECIES_NAME}/2_rawphylo.sh
    bash species/${SPECIES_NAME}/3_multi.sh

done

# clean
find species -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find species -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr

```

