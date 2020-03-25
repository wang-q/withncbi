# Aligning various genera from Protists

Less detailed than *Trichoderma* in
[README.md](https://github.com/wang-q/withncbi/blob/master/pop/README.md).


[TOC levels=1-3]: # ""

- [Aligning various genera from Protists](#aligning-various-genera-from-protists)
  - [Strain info](#strain-info)
  - [NCBI Assembly](#ncbi-assembly)
  - [Count strains](#count-strains)
  - [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
  - [Groups and targets](#groups-and-targets)
  - [Protists: prepare](#protists-prepare)
  - [plasmodium: run](#plasmodium-run)
  - [Protists: run](#protists-run)


## Strain info

| Group          | Genus           | Genus ID | Comments          | Species | Strains |
|:---------------|:----------------|---------:|:------------------|--------:|--------:|
| Alveolata      |                 |          | 囊泡虫类 (顶复虫)    |         |         |
|                | Plasmodium      |     5820 | 疟原虫属           |      20 |      60 |
|                | Toxoplasma      |     5810 | 弓形虫属           |       1 |      16 |
|                | Cryptosporidium |     5806 | 隐孢子虫属          |      13 |      16 |
|                | Eimeria         |     5800 | 艾美球虫           |       2 |       2 |
|                | Theileria       |     5873 | 泰勒虫属           |       4 |       6 |
|                | Babesia         |     5864 | 巴倍虫属           |       6 |       7 |
| Amoebozoa      |                 |          | 变形虫             |         |         |
|                | Acanthamoeba    |     5754 | 棘阿米巴属          |       1 |       2 |
|                | Entamoeba       |     5758 | 内阿米巴属          |       2 |       3 |
|                | Dictyostelium   |     5782 | 网柄菌属           |       4 |       4 |
| Kinetoplastida |                 |          | 动基体目           |         |         |
|                | Leishmania      |     5658 | 利什曼虫属          |      17 |      25 |
|                | Trypanosoma     |     5690 | 锥虫属             |       8 |      13 |
| Stramenopiles  |                 |          | 不等鞭毛类          |         |         |
|                | Blastocystis    |    12967 | 芽囊原虫属          |       2 |       8 |
|                | Nannochloropsis |     5748 | 微拟球藻  (大眼藻纲) |       4 |       8 |
|                | Phytophthora    |     4783 | 疫霉属 (卵菌)       |      11 |      14 |
|                | Pythium         |     4797 | 腐霉属 (卵菌)       |       4 |       4 |
| Euglenozoa     |                 |          | 眼虫门             |         |         |
|                | Crithidia       |     5655 | 短膜虫属           |       2 |       2 |
| Other          |                 |          |                   |         |         |
|                | Giardia         |     5740 | 贾第虫属           |       3 |       6 |

## NCBI Assembly

```bash
export RANK_NAME=Protists

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (
            5820, 5810, 5806, 5800, 5873,
            5864,
            5754, 5758, 5782, 
            5658, 5690,
            12967, 5748, 4783, 4797,
            5655,
            5740
        )
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (
            5820, 5810, 5806, 5800, 5873,
            5864,
            5754, 5758, 5782, 
            5658, 5690,
            12967, 5748, 4783, 4797,
            5655,
            5740
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

# find potential duplicated strains names
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
cd ~/data/alignment/Protists

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Protists.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/Protists.assembly.rsync.sh

bash ASSEMBLY/Protists.assembly.collect.sh

```

## Count strains

```bash
cd ~/data/alignment/Protists

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
        --le 4:3000 \
        --ge 2:100000 |
    tsv-filter -H --ge 3:1000000 |
    tr "\t" "," \
    > ASSEMBLY/n50.pass.csv
        
wc -l ASSEMBLY/n50*
#  204 ASSEMBLY/n50.pass.csv
#  296 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/Protists.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/Protists.assembly.pass.csv

wc -l ASSEMBLY/Protists.assembly*csv
#   293 ASSEMBLY/Protists.assembly.collect.csv
#   201 ASSEMBLY/Protists.assembly.pass.csv

# find potential duplicated strains names
cat ASSEMBLY/Protists.assembly.pass.csv |
    cut -d, -f 7 |
    sort |
    uniq -c |
    sort -nr

```

```bash
cd ~/data/alignment/Protists

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/Protists.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)
    
    n_strains=$(cat ASSEMBLY/Protists.assembly.pass.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        wc -l)
    
    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: $(
        cat ASSEMBLY/Protists.assembly.pass.csv |
            sed -e '1d' |
            cut -d"," -f 2 |
            grep -v "Candidatus" |
            cut -d" " -f 1 |
            sort |
            uniq
    )

#Acanthamoeba    1       2
#Babesia 6       7
#Blastocystis    2       8
#Crithidia       2       2
#Cryptosporidium 13      16
#Dictyostelium   4       4
#Eimeria 2       2
#Entamoeba       2       3
#Giardia 3       6
#Leishmania      17      25
#Nannochloropsis 4       8
#Phytophthora    11      14
#Plasmodium      20      60
#Pythium 4       4
#Theileria       4       6
#Toxoplasma      1       16
#Trypanosoma     8       13

```

## Raw phylogenetic tree by MinHash

```bash
mkdir -p ~/data/alignment/Protists/mash
cd ~/data/alignment/Protists/mash

for name in $(cat ../ASSEMBLY/Protists.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 ); do
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
    cat ../ASSEMBLY/Protists.assembly.pass.csv | sed -e '1d' | cut -d"," -f 1 | parallel echo "{}.msh"
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
    rsvg-convert -o ~/Scripts/withncbi/image/Protists.png

```

## Groups and targets

Review `ASSEMBLY/Protists.assembly.pass.csv` and `mash/groups.tsv`

| #Serial | Group              | Count | Target                 | Sequencing     |
|:--------|:-------------------|:------|:-----------------------|:---------------|
| 1       | A_Ei               | 4     | A_cas                  |                |
| 2       | Babesia            | 7     | Ba_bov_T2Bo            | 8X Sanger, WGS |
| 4       | Blastocystis       | 8     | Bl_hom                 |                |
| 8       | Cryptosporidium    | 11    | Cry_parvum_Iowa_II     |                |
| 7       | Dictyostelium      | 11    | D_disc_AX4             |                |
| 10      | Giardia            | 6     | G_intes                |                |
| 11      | Leishmania         | 24    | L_maj_Friedlin         |                |
| 12      | Nannochloropsis    | 8     | N_oce                  |                |
| 13      | Phytophthora       | 10    | Ph_soj                 |                |
| 15      | Pl_ber_cha_vin_yoe | 10    | Pl_yoe                 |                |
| 16      | Pl_kno_viv         | 12    | Pl_viv                 |                |
| 17      | Pl_falcip          | 31    | Pl_falcip_3D7          |                |
| 18      | Pythium            | 4     | Py_gui                 |                |
| 19      | Theileria          | 6     | Th_parva_Muguga        |                |
| 20      | Toxoplasma         | 16    | To_gondii_ME49         |                |
| 21      | Trypanosoma        | 5     | Tr_bruc_brucei_TREU927 |                |

```bash
mkdir -p ~/data/alignment/Protists/taxon
cd ~/data/alignment/Protists/taxon

cp ../mash/tree.nwk .

# manually combine Ba_mic
# manually remove bad assemblies
cat ../mash/groups.tsv |
    grep -v "Ba_mic" |
    grep -v "Cry_bai" |
    grep -v "En_inv_IP1" |
    grep -v "L_sp_A" \
    > groups.tsv
echo -e "2\tBa_mic" >> groups.tsv
echo -e "2\tBa_mic_RI" >> groups.tsv
#echo -e "11\tCr_win_CBS_7118" >> groups.tsv
#echo -e "11\tCr_dep_CBS_7841" >> groups.tsv
#echo -e "11\tCr_dep_CBS_7855" >> groups.tsv

ARRAY=(
    'A_Ei::A_cas'
    'Babesia::Ba_bov_T2Bo' 
    'Blastocystis::Bl_hom' 
    'Cryptosporidium::Cry_parvum_Iowa_II' 
    'Dictyostelium::D_disc_AX4' 
    'Giardia::G_intes' 
    'Leishmania::L_maj_Friedlin' 
    'Nannochloropsis::N_oce' 
    'Phytophthora::Ph_soj' 
    'Pl_ber_cha_vin_yoe::Pl_yoe' 
    'Pl_kno_viv::Pl_viv' 
    'Pl_falcip::Pl_falcip_3D7' 
    'Pythium::Py_gui' 
    'Theileria::Th_parva_Muguga'
    'Toxoplasma::To_gondii_ME49'
    'Trypanosoma::Tr_bruc_brucei_TREU927'
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
Ba_big
Ba_bov_T2Bo
Ba_mic_RI
Cry_parvum_Iowa_II
D_disc_AX4
L_braz_MHOM_BR_75_M2904
L_don
L_infa_JPCM5
L_maj_Friedlin
L_mex_MHOM_GT_2001_U1103
Pl_falcip_3D7
Pl_kno_H
Pl_viv
Th_ann
Th_equi_WA
Th_ori_Shintoku
Th_parva
Th_parva_Muguga
Tr_bruc_brucei_TREU927
Tr_bruc_gambiense_DAL972
EOF

```

## Protists: prepare

* Rsync to hpcc

```bash
rsync -avP \
    ~/data/alignment/Protists/ \
    wangq@202.119.37.251:data/alignment/Protists

# rsync -avP wangq@202.119.37.251:data/alignment/Protists/ ~/data/alignment/Protists

```

* `--perseq` for Chromosome-level assemblies and targets
* Use `Stramenopiles` as `Eukaryota` takes too long to compute

```bash
cd ~/data/alignment/Protists/

$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Eukaryota -stat
#174 ancestral and ubiquitous sequence(s) with a total length of 45804 bp
#0 eukaryota specific repeats with a total length of 0 bp
#31314 lineage specific sequence(s) with a total length of 81842705 bp

$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Alveolata -stat
#174 ancestral and ubiquitous sequence(s) with a total length of 45804 bp
#3 alveolata specific repeats with a total length of 7681 bp
#30 lineage specific sequence(s) with a total length of 81127 bp

$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Kinetoplastida -stat
#176 ancestral and ubiquitous sequence(s) with a total length of 48931 bp
#0 kinetoplastida specific repeats with a total length of 0 bp
#44 lineage specific sequence(s) with a total length of 67956 bp

$(brew --prefix repeatmasker)/libexec/util/queryRepeatDatabase.pl \
    -species Stramenopiles -stat
#174 ancestral and ubiquitous sequence(s) with a total length of 45804 bp
#3 stramenopiles specific repeats with a total length of 6250 bp
#620 lineage specific sequence(s) with a total length of 2090320 bp

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    $( cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 | parallel -j 1 echo " --perseq {} " ) \
    $( cat taxon/chr-level.list | parallel -j 1 echo " --perseq {} " ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Stramenopiles --parallel 24"

bsub -q mpi -n 24 -J "Protists-0_prep" "bash GENOMES/0_prep.sh"

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


## plasmodium: run

Plasmodium distributed on many branches.

```bash
cd ~/data/alignment/Protists/

# sanger
egaz template \
    GENOMES/Pl_falcip_3D7 \
    GENOMES/Pl_yoe \
    GENOMES/Pl_viv \
    GENOMES/Pl_ber_ANKA \
    GENOMES/Pl_cha_chabaudi \
    GENOMES/Pl_cyn_B \
    GENOMES/Pl_kno_H \
    GENOMES/Pl_ovale \
    GENOMES/Pl_rel \
    --multi -o groups/plasmodium/ \
    --multiname sanger \
    --tree taxon/tree.nwk \
    --parallel 24 -v

bsub -q mpi -n 24 -J "plasmodium-1_pair" "bash groups/plasmodium/1_pair.sh"
bsub  -w "ended(plasmodium-1_pair)" \
    -q mpi -n 24 -J "plasmodium-3_multi" "bash groups/plasmodium/3_multi.sh"

```

## Protists: run

```bash
cd ~/data/alignment/Protists/

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

