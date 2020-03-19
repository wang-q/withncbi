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
  - [Ca_alb: run](#ca_alb-run)

## Strain info

| Group          | Genus            | Genus ID | Comments          | Species | Strains |
|:---------------|:-----------------|---------:|:------------------|--------:|--------:|
| Ascomycetes    |                  |          |                   |         |         |
|                | Saccharomyces    |     4930 | 酵母菌属           |       8 |      50 |
|                | Aspergillus      |     5052 | 曲霉菌属           |      47 |      62 |
|                | Blastomyces      |   229219 | 芽生菌属 (lung)    |       2 |       4 |
|                | Candida          |  1535326 | 念珠菌属           |       4 |      37 |
|                | Coccidioides     |     5500 | 球孢子菌属 (lung)   |       2 |      14 |
|                | Colletotrichum   |     5455 | 炭疽菌属           |       6 |       6 |
|                | Epichloe         |     5112 |                   |       2 |       3 |
|                | Fusarium         |     5506 | 镰刀菌             |       6 |      40 |
|                | Hanseniaspora    |    29832 | 有孢汉逊酵母        |       3 |       4 |
|                | Histoplasma      |     5036 | 组织胞浆菌属 (lung) |       1 |       5 |
|                | Kazachstania     |    71245 |                   |       2 |       2 |
|                | Metschnikowia    |    27320 | 梅奇酵母属          |       3 |       3 |
|                | Ogataea          |   461281 |                   |       1 |       1 |
|                | Paracoccidioides |    38946 | 副球孢子菌属 (lung) |       2 |       3 |
|                | Penicillium      |     5073 | 青霉菌属           |      12 |      18 |
|                | Pichia           |     4919 | 毕赤酵母属          |       2 |       2 |
|                | Pneumocystis     |     4753 | 肺孢子菌属 (lung)   |       3 |       3 |
|                | Pyricularia      |    48558 | 梨孢属             |       2 |       7 |
|                | Sporothrix       |    29907 | 孢子丝菌属 (skin)   |       3 |       4 |
|                | Talaromyces      |     5094 | 踝节菌属 (lung)    |       2 |       3 |
|                | Trichoderma      |     5543 | 木霉属             |       7 |      10 |
|                | Trichophyton     |     5550 | 毛癣菌属           |       7 |      16 |
|                | Verticillium     |  1036719 | 轮枝菌属           |       3 |       6 |
|                | Yarrowia         |     4951 | 耶氏酵母           |       1 |       3 |
|                | Zymoseptoria     |  1047167 |                   |       4 |      16 |
| Basidiomycetes |                  |          |                   |         |         |
|                | Cryptococcus     |     5206 | 隐球菌属, 脑膜炎    |       5 |      53 |
|                | Malassezia       |    55193 | 马拉色菌属          |       3 |       3 |
|                | Puccinia         |     5296 | 柄锈菌属           |       4 |       6 |
|                | Rhodotorula      |     5533 | 红酵母属           |       2 |       4 |
|                | Ustilago         |     5269 | 黑粉菌属           |       1 |       1 |
| Basidiomycetes |                  |          |                   |         |         |
|                | Mucor            |     4830 | 毛霉菌属           |       2 |       4 |

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

mkdir -p taxon

parallel --no-run-if-empty --linebuffer -k -j 4 '
    cat ASSEMBLY/Fungi.assembly.pass.csv |
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
    > genus.list

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
        mash sketch -k 21 -s 100000 -p 4 - -I "${name}" -o ${name}
done

mash triangle -E -p 4 -l <( find . -maxdepth 1 -type f -name "*.msh" | sort ) > dist.tsv

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

nw_display -w 600 -b 'visibility:hidden' -s tree.nwk |
    rsvg-convert -o ~/Scripts/withncbi/image/Fungi.png

```

## Groups and targets

Review `ASSEMBLY/Fungi.assembly.pass.csv` and `mash/groups.tsv`

| #Serial | Group             | Count | Target                                   | Sequencing technology  |
|:--------|:------------------|------:|:-----------------------------------------|:-----------------------|
| 4       | A_fum             |     9 | A_fum_Af293                              | 10.5X Sanger, WGS      |
| 2       | A_nig             |    22 | A_nig_CBS_513_88                         | 7.5X Sanger, BAC       |
| 5       | A_ory             |    10 | A_ory_RIB40                              | 9X Sanger, WGS         |
| 3       | A_nid             |    12 | A_nid_FGSC_A4                            | 13X Sanger, BAC        |
| 6       | B_Hi_Pa           |    12 | B_der_ER_3                               | 9.4X Sanger, WGS       |
|         |                   |       | Hi_cap_G186AR                            | 10X Sanger             |
|         |                   |       | Pa_bras_Pb18                             | 9.8X Sanger, WGS       |
| 7       | Ca_alb            |    35 | Ca_alb_WO_1                              | 10X Sanger, WGS        |
| 8       | Ca_Ha_K_Me_O_Pi_U |    15 | U_may_521                                | 10X Sanger, WGS        |
|         |                   |       | Ca_ort_Co_90_125                         | 10X 454, fosmid        |
|         |                   |       | K_afr_CBS_2517                           | 20X 454                |
| 9       | Coccidioides      |    14 | Coc_imm_RS                               | 14.4X Sanger, WGS      |
| 10      | Colletotrichum    |     6 | Col_hig_IMI_349063                       | 133.0X PacBio          |
| 11      | Cr_gattii         |    17 | Cr_gattii_WM276                          | 6X Sanger, fosmid, BAC |
| 14      | Cr_neof           |    31 | Cr_neof_var_grubii_H99                   | 11X Sanger, WGS        |
|         |                   |       | Cr_neof_var_neoformans_JEC21             | 12.5X Sanger, WGS      |
| 15      | Epichloe          |     3 | E_fes_Fl1                                | 210.0X PacBio          |
| 16      | Fusarium          |    40 | F_oxy_f_sp_lycopersici_4287              | 6X Sanger, WGS         |
|         |                   |       | F_gramine_PH_1                           | 10X Sanger, WGS        |
|         |                   |       | F_vert_7600                              | 8X Sanger, WGS         |
| 17      | Malassezia        |     3 | Ma_glob_CBS_7966                         | 7X Sanger              |
| 18      | Mucor             |     4 | Mu_lus_CBS_277_49                        | 9.49X Sanger           |
| 19      | Penicillium       |    14 | Pe_rubens_Wisconsin_54_1255              | 9.8X Sanger, WGS       |
| 21      | Pneumocystis      |     3 | Pn_jir_RU7                               | 241.4X Illumina        |
| 22      | Puccinia          |     6 | Pu_graminis_f_sp_tritici_CRL_75_36_700_3 | 7.88X Sanger, WGS      |
| 23      | Pyricularia       |     7 | Py_ory_70_15                             | 7X Sanger, WGS         |
| 24      | Rhodotorula       |     4 | R_graminis_WP1                           | 8.55X Sanger, WGS      |
| 25      | Saccharomyces     |    50 | Sa_cerevisiae_S288C                      | Reference              |
| 26      | Sporothrix        |     4 | Sp_ins_RCEF_264                          | 7037.0X HiSeq          |
| 27      | Talaromyces       |     3 | Ta_mar_ATCC_18224                        | 8.8X Sanger            |
| 28      | Trichoderma       |    10 | Trichod_atr_IMI_206040                   | 8.26X Sanger           |
| 29      | Trichophyton      |    16 | Trichop_rubr_CBS_118892                  | 8.19X Sanger, WGS      |
| 30      | Verticillium      |     6 | V_alf_VaMs_102                           | 4.08X Sanger, WGS      |
| 31      | Yarrowia          |     3 | Y_lip_CLIB122                            | 10X Sanger, WGS        |
| 32      | Zymoseptoria      |    16 | Z_tritici_IPO323                         | 8.9X Sanger, WGS       |

```bash
mkdir -p ~/data/alignment/Fungi/taxon
cd ~/data/alignment/Fungi/taxon

cp ../mash/tree.nwk .
cp ../mash/groups.tsv .

ARRAY=(
    'A_fum::A_fum_Af293' # group 4
    'A_nig::A_nig_CBS_513_88' # group 2
    'A_ory::A_ory_RIB40' # group 5
    'A_nid::A_nid_FGSC_A4' # group 3
    'B_Hi_Pa::B_der_ER_3' # 6
    'Ca_alb::Ca_alb_WO_1' # 7
    'Ca_Ha_K_Me_O_Pi_U::U_may_521' # 8
    'Coccidioides::Coc_imm_RS' # 9
    'Colletotrichum::Col_hig_IMI_349063' # 10
    'Cr_gattii::Cr_gattii_WM276' # 11
    'Cr_neof::Cr_neof_var_grubii_H99' # 14
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

mlr --itsv --omd cat group_target.tsv

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
    --perseq A_fum_Af293 \
    --perseq A_ory_RIB40 \
    --perseq Ca_alb_WO_1 \
    --perseq Ca_dub_CD36 \
    --perseq Ca_ort_Co_90_125 \
    --perseq Cr_gattii_VGI_Cryptococcus_gattii_WM276 \
    --perseq Cr_neof_var_grubii_H99 \
    --perseq Cr_neof_var_neoformans_B_3501A \
    --perseq Cr_neof_var_neoformans_JEC21 \
    --perseq F_fuj_IMI_58289 \
    --perseq F_gramine_PH_1 \
    --perseq F_oxy_f_sp_lycopersici_4287 \
    --perseq F_vert_7600 \
    --perseq K_afr_CBS_2517 \
    --perseq K_nag_CBS_8797 \
    --perseq Py_ory_70_15 \
    --perseq Sa_cerevisiae_S288C \
    --perseq Sa_cerevisiae_Sigma1278b \
    --perseq U_may_521 \
    --perseq Y_lip_CLIB122 \
    --perseq Y_lip_WSH_Z06 \
    --perseq Z_tritici_IPO323 \
    --perseq Z_tritici_ST99CH_1A5 \
    --perseq Z_tritici_ST99CH_1E4 \
    --perseq Z_tritici_ST99CH_3D1 \
    --perseq Z_tritici_ST99CH_3D7 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Fungi --parallel 24"

bsub -q mpi -n 24 -J "Fungi-0_prep" "bash GENOMES/0_prep.sh"

ls -t output.* | head -n 1 | xargs tail -f | grep "==>"

# gff
for n in Calb_WO_1 Cdub_CD36 Ctro_MYA_3404; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"
    
    gzip -d -c ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

## Ca_alb: run

```bash
# sanger
egaz template \
    GENOMES/Calb_WO_1 \
    GENOMES/Cdub_CD36 \
    GENOMES/Ctro_MYA_3404 \
    --multi -o multi/ \
    --multiname sanger --order \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_pair" "bash multi/1_pair.sh"
bsub  -w "ended(candida-1_pair)" \
    -q mpi -n 24 -J "candida-3_multi" "bash multi/3_multi.sh"

# multi
egaz template \
    GENOMES/Calb_WO_1 \
    $(find GENOMES -maxdepth 1 -mindepth 1 -type d | grep -v "Calb_WO_1") \
    --multi -o multi/ \
    --tree mash/tree.nwk \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_pair" "bash multi/1_pair.sh"
bsub  -w "ended(candida-1_pair)" \
    -q mpi -n 24 -J "candida-3_multi" "bash multi/3_multi.sh"

# multi_Calb
egaz template \
    GENOMES/Calb_WO_1 \
    $(find GENOMES -maxdepth 1 -mindepth 1 -type d | grep "Calb_" | grep -v "Calb_SC5314") \
    GENOMES/Cdub_CD36 \
    --multi -o multi/ \
    --multiname Calb --tree mash/tree.nwk --outgroup Cdub_CD36 \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-3_multi" "bash multi/3_multi.sh"

# self
egaz template \
    GENOMES/Calb_WO_1 \
    GENOMES/Cdub_CD36 \
    GENOMES/Ctro_MYA_3404 \
    --self -o self/ \
    --circos --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_self" "bash self/1_self.sh"
bsub -w "ended(candida-1_self)" \
    -q mpi -n 24 -J "candida-3_proc" "bash self/3_proc.sh"
bsub  -w "ended(candida-3_proc)" \
    -q mpi -n 24 -J "candida-4_circos" "bash self/4_circos.sh"

```

