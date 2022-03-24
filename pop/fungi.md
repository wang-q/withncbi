# Aligning various genera from Fungi

Less detailed than *Trichoderma* in
[README.md](https://github.com/wang-q/withncbi/blob/master/pop/README.md).


[TOC levels=1-3]: # ""

- [Aligning various genera from Fungi](#aligning-various-genera-from-fungi)
    - [Strain info](#strain-info)
    - [NCBI Assembly](#ncbi-assembly)
    - [Count strains](#count-strains)
    - [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
    - [Groups and targets](#groups-and-targets)
    - [Fungi: prepare](#fungi-prepare)
    - [candida: run](#candida-run)
    - [Fungi: run](#fungi-run)

## Strain info

| Group          | Genus            | Genus ID | Species | Strains | Comments      |
|:---------------|:-----------------|---------:|--------:|--------:|:--------------|
| Ascomycetes    |                  |          |         |         |               |
|                | Saccharomyces    |     4930 |       8 |      50 | 酵母菌属          |
|                | Aspergillus      |     5052 |      47 |      62 | 曲霉菌属          |
|                | Blastomyces      |   229219 |       2 |       4 | 芽生菌属 (lung)   |
|                | Candida          |  1535326 |       4 |      37 | 念珠菌属          |
|                | Coccidioides     |     5500 |       2 |      14 | 球孢子菌属 (lung)  |
|                | Colletotrichum   |     5455 |       6 |       6 | 炭疽菌属          |
|                | Epichloe         |     5112 |       2 |       3 |               |
|                | Fusarium         |     5506 |       6 |      40 | 镰刀菌           |
|                | Hanseniaspora    |    29832 |       3 |       4 | 有孢汉逊酵母        |
|                | Histoplasma      |     5036 |       1 |       5 | 组织胞浆菌属 (lung) |
|                | Kazachstania     |    71245 |       2 |       2 |               |
|                | Metschnikowia    |    27320 |       3 |       3 | 梅奇酵母属         |
|                | Ogataea          |   461281 |       1 |       1 |               |
|                | Paracoccidioides |    38946 |       2 |       3 | 副球孢子菌属 (lung) |
|                | Penicillium      |     5073 |      12 |      18 | 青霉菌属          |
|                | Pichia           |     4919 |       2 |       2 | 毕赤酵母属         |
|                | Pneumocystis     |     4753 |       3 |       3 | 肺孢子菌属 (lung)  |
|                | Pyricularia      |    48558 |       2 |       7 | 梨孢属           |
|                | Sporothrix       |    29907 |       3 |       4 | 孢子丝菌属 (skin)  |
|                | Talaromyces      |     5094 |       2 |       3 | 踝节菌属 (lung)   |
|                | Trichoderma      |     5543 |       7 |      10 | 木霉属           |
|                | Trichophyton     |     5550 |       7 |      16 | 毛癣菌属          |
|                | Verticillium     |  1036719 |       3 |       6 | 轮枝菌属          |
|                | Yarrowia         |     4951 |       1 |       3 | 耶氏酵母          |
|                | Zymoseptoria     |  1047167 |       4 |      16 |               |
| Basidiomycetes |                  |          |         |         |               |
|                | Cryptococcus     |     5206 |       5 |      53 | 隐球菌属, 脑膜炎     |
|                | Malassezia       |    55193 |       3 |       3 | 马拉色菌属         |
|                | Puccinia         |     5296 |       4 |       6 | 柄锈菌属          |
|                | Rhodotorula      |     5533 |       2 |       4 | 红酵母属          |
|                | Ustilago         |     5269 |       1 |       1 | 黑粉菌属          |
| Other          |                  |          |         |         |               |
|                | Mucor            |     4830 |       2 |       4 | 毛霉菌属          |

* https://patient.info/doctor/fungal-lung-infections

* https://www.cdc.gov/fungal/diseases/index.html

## NCBI Assembly

```bash
export RANK_NAME=Fungi

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND genus_id in (
            4930, 5052, 229219, 1535326, 5500,
            5455, 5112, 5506, 29832, 5036,
            71245, 27320, 461281, 38946, 5073,
            4753, 4919, 48558, 29907, 5094,
            5543, 5550, 1036719, 4951, 1047167,
            5206, 5296, 55193, 5533, 5269,
            4830
        )
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND genus_id in (
            4930, 5052, 229219, 1535326, 5500,
            5455, 5112, 5506, 29832, 5036,
            71245, 27320, 461281, 38946, 5073,
            4753, 4919, 48558, 29907, 5094,
            5543, 5550, 1036719, 4951, 1047167,
            5206, 5296, 55193, 5533, 5269,
            4830
        )
    " \
    >> raw.tsv

cat raw.tsv |
    grep -v '^#' |
    tsv-filter --not-regex "1:\[.+\]" |
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
cd ~/data/alignment/Fungi

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Fungi.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/Fungi.assembly.rsync.sh

bash ASSEMBLY/Fungi.assembly.collect.sh

```

## Count strains

```bash
cd ~/data/alignment/Fungi

for dir in $(find ASSEMBLY -maxdepth 1 -mindepth 1 -type d | sort); do
    1>&2 echo "==> ${dir}"
    name=$(basename ${dir})

    find ${dir} -type f -name "*_genomic.fna.gz" |
        grep -v "_from_" |
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
        --le 4:2000 \
        --ge 2:100000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv

wc -l ASSEMBLY/n50*
#  401 ASSEMBLY/n50.pass.csv
#  453 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/Fungi.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/Fungi.assembly.pass.csv

wc -l ASSEMBLY/Fungi.assembly*csv
#   445 ASSEMBLY/Fungi.assembly.collect.csv
#   394 ASSEMBLY/Fungi.assembly.pass.csv

```

```bash
cd ~/data/alignment/Fungi

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/Fungi.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)

    n_strains=$(cat ASSEMBLY/Fungi.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        wc -l)

    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: $(
        cat ASSEMBLY/Fungi.assembly.pass.csv |
            sed -e '1d' |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            cut -d" " -f 1 |
            sort |
            uniq
    )

#Aspergillus     47      62
#Blastomyces     2       4
#Candida 4       37
#Coccidioides    2       14
#Colletotrichum  6       6
#Cryptococcus    5       53
#Epichloe        2       3
#Fusarium        6       40
#Hanseniaspora   3       4
#Histoplasma     1       5
#Kazachstania    2       2
#Malassezia      3       3
#Metschnikowia   3       3
#Mucor   2       4
#Ogataea 1       1
#Paracoccidioides        2       3
#Penicillium     12      18
#Pichia  2       2
#Pneumocystis    3       3
#Puccinia        4       6
#Pyricularia     2       7
#Rhodotorula     2       4
#Saccharomyces   8       50
#Sporothrix      3       4
#Talaromyces     2       3
#Trichoderma     7       10
#Trichophyton    7       16
#Ustilago        1       1
#Verticillium    3       6
#Yarrowia        1       3
#Zymoseptoria    4       16


```

## Raw phylogenetic tree by MinHash

```bash
mkdir -p ~/data/alignment/Fungi/mash
cd ~/data/alignment/Fungi/mash

for name in $(cat ../ASSEMBLY/Fungi.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 ); do
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
    cat ../ASSEMBLY/Fungi.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 | parallel echo "{}.msh"
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
    rsvg-convert -o ~/Scripts/withncbi/image/Fungi.png

```

## Groups and targets

Review `ASSEMBLY/Fungi.assembly.pass.csv` and `mash/groups.tsv`

| #Serial | Group             | Count | Target                      | Sequencing             |
|:--------|:------------------|------:|:----------------------------|:-----------------------|
| 1       | A_aculeati        |    10 | A_aculeati_CBS_121060       | 95.8X Illumina         |
| 2       | A_nig             |    22 | A_nig_CBS_513_88            | 7.5X Sanger, BAC       |
| 3       | A_nid             |    12 | A_nid_FGSC_A4               | 13X Sanger, BAC        |
| 4       | A_fum             |     9 | A_fum_Af293                 | 10.5X Sanger, WGS      |
| 5       | A_ory             |    10 | A_ory_RIB40                 | 9X Sanger, WGS         |
| 6       | B_Hi_Pa           |    12 | B_der_ER_3                  | 9.4X Sanger, WGS       |
|         |                   |       | Hi_cap_G186AR               | 10X Sanger             |
|         |                   |       | Pa_bras_Pb18                | 9.8X Sanger, WGS       |
| 7       | Ca_alb            |    35 | Ca_alb_WO_1                 | 10X Sanger, WGS        |
| 8       | Ca_Ha_K_Me_O_Pi_U |    15 | U_may_521                   | 10X Sanger, WGS        |
|         |                   |       | Ca_ort_Co_90_125            | 10X 454, fosmid        |
|         |                   |       | K_afr_CBS_2517              | 20X 454                |
| 9       | Coccidioides      |    14 | Coc_imm_RS                  | 14.4X Sanger, WGS      |
| 10      | Colletotrichum    |     6 | Col_hig_IMI_349063          | 133.0X PacBio          |
| 11      | Cr_amy            |     5 | Cr_amy_CBS_6039             | 117.0X Illumina        |
| 12      | Cr_gattii         |    17 | Cr_gattii_WM276             | 6X Sanger, fosmid, BAC |
| 14      | Cr_neof           |    31 | Cr_neof_grubii_H99          | 11X Sanger, WGS        |
|         |                   |       | Cr_neof_neoformans_JEC21    | 12.5X Sanger, WGS      |
| 15      | Epichloe          |     3 | E_fes_Fl1                   | 210.0X PacBio          |
| 16      | Fusarium          |    40 | F_oxy_f_sp_lycopersici_4287 | 6X Sanger, WGS         |
|         |                   |       | F_gramine_PH_1              | 10X Sanger, WGS        |
|         |                   |       | F_vert_7600                 | 8X Sanger, WGS         |
| 17      | Malassezia        |     3 | Ma_glob_CBS_7966            | 7X Sanger              |
| 18      | Mucor             |     4 | Mu_lus_CBS_277_49           | 9.49X Sanger           |
| 19      | Penicillium       |    14 | Pe_rubens_Wisconsin_54_1255 | 9.8X Sanger, WGS       |
| 21      | Pneumocystis      |     3 | Pn_jir_RU7                  | 241.4X Illumina        |
| 22      | Puccinia          |     6 | Pu_graminis_CRL_75_36_700_3 | 7.88X Sanger, WGS      |
| 23      | Pyricularia       |     7 | Py_ory_70_15                | 7X Sanger, WGS         |
| 24      | Rhodotorula       |     4 | R_graminis_WP1              | 8.55X Sanger, WGS      |
| 25      | Saccharomyces     |    50 | Sa_cerevisiae_S288C         | Reference              |
| 26      | Sporothrix        |     4 | Sp_ins_RCEF_264             | 7037.0X HiSeq          |
| 27      | Talaromyces       |     3 | Ta_mar_ATCC_18224           | 8.8X Sanger            |
| 28      | Trichoderma       |    10 | Trichod_atr_IMI_206040      | 8.26X Sanger           |
| 29      | Trichophyton      |    16 | Trichop_rubr_CBS_118892     | 8.19X Sanger, WGS      |
| 30      | Verticillium      |     6 | V_alf_VaMs_102              | 4.08X Sanger, WGS      |
| 31      | Yarrowia          |     3 | Y_lip_CLIB122               | 10X Sanger, WGS        |
| 32      | Zymoseptoria      |    16 | Z_tritici_IPO323            | 8.9X Sanger, WGS       |

```bash
mkdir -p ~/data/alignment/Fungi/taxon
cd ~/data/alignment/Fungi/taxon

cp ../mash/tree.nwk .

# manually combine Cr_amy and Cr_dep
cat ../mash/groups.tsv |
    grep -v "Cr_amy_" |
    grep -v "Cr_win_" |
    grep -v "Cr_dep_" \
    > groups.tsv
echo -e "11\tCr_amy_CBS_6039" >> groups.tsv
echo -e "11\tCr_amy_CBS_6273" >> groups.tsv
echo -e "11\tCr_win_CBS_7118" >> groups.tsv
echo -e "11\tCr_dep_CBS_7841" >> groups.tsv
echo -e "11\tCr_dep_CBS_7855" >> groups.tsv

ARRAY=(
    'A_aculeati::A_aculeati_CBS_121060' # group 1
    'A_nig::A_nig_CBS_513_88' # group 2
    'A_nid::A_nid_FGSC_A4' # group 3
    'A_fum::A_fum_Af293' # group 4
    'A_ory::A_ory_RIB40' # group 5
    'B_Hi_Pa::B_der_ER_3' # 6
    'Ca_alb::Ca_alb_WO_1' # 7
    'Ca_Ha_K_Me_O_Pi_U::U_may_521' # 8
    'Coccidioides::Coc_imm_RS' # 9
    'Colletotrichum::Col_hig_IMI_349063' # 10
    'Cr_amy::Cr_amy_CBS_6039' # 11
    'Cr_gattii::Cr_gattii_WM276' # 12
    'Cr_neof::Cr_neof_grubii_H99' # 14
    'Epichloe::E_fes_Fl1' # 15
    'Fusarium::F_oxy_f_sp_lycopersici_4287' # 16
    'Malassezia::Ma_glob_CBS_7966' # 17
    'Mucor::Mu_lus_CBS_277_49' # 18
    'Penicillium::Pe_rubens_Wisconsin_54_1255' # 19
    'Pneumocystis::Pn_jir_RU7' # 21
    'Puccinia::Pu_graminis_CRL_75_36_700_3' # 22
    'Pyricularia::Py_ory_70_15' # 23
    'Rhodotorula::R_graminis_WP1' # 24
    'Saccharomyces::Sa_cerevisiae_S288C' # 25
    'Sporothrix::Sp_ins_RCEF_264' # 26
    'Talaromyces::Ta_mar_ATCC_18224' # 27
    'Trichoderma::Trichod_atr_IMI_206040' # 28
    'Trichophyton::Trichop_rubr_CBS_118892' # 29
    'Verticillium::V_alf_VaMs_102' # 30
    'Yarrowia::Y_lip_CLIB122' # 31
    'Zymoseptoria::Z_tritici_IPO323' # 32
)

echo -e "#Serial\tGroup\tCount\tTarget\tSequencing" > group_target.tsv

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

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}\t" >> group_target.tsv

done

mlr --itsv --omd cat group_target.tsv

cat <<'EOF' > chr-level.list
Ca_dub_CD36
Ca_ort_Co_90_125
Cr_neof_neoformans_JEC21
F_fuj_IMI_58289
F_gramine_PH_1
F_vert_7600
Hi_cap_G186AR
K_afr_CBS_2517
Pa_bras_Pb18
Sa_cerevisiae_Sigma1278b
EOF

```

## Fungi: prepare

* Rsync to hpcc

```bash
rsync -avP \
    ~/data/alignment/Fungi/ \
    wangq@202.119.37.251:data/alignment/Fungi

# rsync -avP wangq@202.119.37.251:data/alignment/Fungi/ ~/data/alignment/Fungi

```

`--perseq` for Chromosome-level assemblies and targets

```bash
cd ~/data/alignment/Fungi/

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    $( cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 | parallel -j 1 echo " --perseq {} " ) \
    $( cat taxon/chr-level.list | parallel -j 1 echo " --perseq {} " ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Fungi --parallel 24"

bsub -q mpi -n 24 -J "Fungi-0_prep" "bash GENOMES/0_prep.sh"

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

## candida: run

`Ca_ort_Co_90_125` is on another branch

```bash
cd ~/data/alignment/Fungi/

# sanger
egaz template \
    GENOMES/Ca_alb_WO_1 \
    GENOMES/Ca_dub_CD36 \
    GENOMES/Ca_tro_MYA_3404 \
    GENOMES/Ca_ort_Co_90_125 \
    --multi -o groups/candida/ \
    --multiname sanger \
    --tree taxon/tree.nwk \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_pair" "bash groups/candida/1_pair.sh"
bsub  -w "ended(candida-1_pair)" \
    -q mpi -n 24 -J "candida-3_multi" "bash groups/candida/3_multi.sh"

# multi
egaz template \
    GENOMES/Ca_alb_WO_1 \
    $(cat taxon/Ca_alb | grep -v -x "Ca_alb_WO_1" | parallel -j 1 echo "GENOMES/{}") \
    --multi -o groups/candida/ \
    --tree taxon/tree.nwk \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_pair" "bash groups/candida/1_pair.sh"
bsub -w "ended(candida-1_pair)" \
    -q mpi -n 24 -J "candida-3_multi" "bash groups/candida/3_multi.sh"

```

## Fungi: run

```bash
cd ~/data/alignment/Fungi/

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

```

