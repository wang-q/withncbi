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

| Group                   | Species               | Species ID | Comments      | Strains |
|:------------------------|:----------------------|-----------:|:--------------|--------:|
| M. tuberculosis complex |                       |            |               |         |
|                         | M. canetti            |      78331 | 卡氏分枝杆菌    |       5 |
|                         | M. decipiens          |    1430326 |               |         |
|                         | M. mungi              |    1844474 |               |         |
|                         | M. orygis             |    1305738 |               |       0 |
|                         | M. tuberculosis       |       1773 | 结核分枝杆菌    |     187 |
| M. avium complex        |                       |            |               |         |
|                         | M. avium              |       1764 | 鸟分枝杆菌      |       7 |
|                         | M. chimaera           |     222805 |               |         |
|                         | M. colombiense        |     339268 |               |       1 |
|                         | M. intracellulare     |       1767 | 胞内分枝杆菌    |       8 |
|                         | M. mantenii           |     560555 |               |         |
|                         | M. marseillense       |     701042 |               |         |
|                         | M. paraintracellulare |    1138383 | 副胞内分枝杆菌   |         |
| M. simiae complex       |                       |            |               |         |
|                         | M. genavense          |      36812 |               |       1 |
|                         | M. parascrofulaceum   |     240125 |               |       1 |
|                         | M. simiae             |       1784 | 猿分支杆菌      |       1 |
| Close to MSC            |                       |            |               |         |
|                         | M. asiaticum          |       1790 | 亚洲分枝杆菌    |       1 |
|                         | M. bohemicum          |       1998 | 波希米亚分枝杆菌 |       1 |
| Others                  |                       |            |               |         |
|                         | M. gordonae           |       1778 | 戈登分枝杆菌    |         |
|                         | M. haemophilum        |      29311 |               |       1 |
|                         | M. kansasii           |       1955 | 堪萨斯分枝杆菌   |       4 |
|                         | M. leprae             |       1769 | 麻风分枝杆菌    |       3 |
|                         | M. liflandii          |     261524 |               |       1 |
|                         | M. malmoense          |       1780 | 莫尔门分枝杆菌   |         |
|                         | M. marinum            |       1781 | 海洋分枝杆菌    |       4 |
|                         | M. persicum           |    1487726 |               |         |
|                         | M. pseudoshottsii     |     265949 |               |       1 |
|                         | M. ulcerans           |       1809 |               |       2 |
|                         | M. xenopi             |       1789 | 蟾分枝杆菌      |       3 |

| Subspecies                      |     ID | Comments        | Strains |
|:--------------------------------|-------:|:----------------|--------:|
| M. tuberculosis var. africanum  |  33894 | 非洲分枝杆菌      |      26 |
| M. tuberculosis var. bovis      |   1765 | 牛分枝杆菌, 卡介苗 |      14 |
| M. tuberculosis var. caprae     | 115862 |                 |       1 |
| M. tuberculosis var. microti    |   1806 | 仓鼠分枝杆菌      |       1 |
| M. tuberculosis var. pinnipedii | 194542 |                 |       1 |

M. tuberculosis var. microti is distinct from other M. tuberculosis strains

* Other species of Actinobacteria

| Species                     | Species ID | Comments       | Strains |
|:----------------------------|-----------:|:---------------|--------:|
| Amycolatopsis mediterranei  |      33910 | 地中海拟无枝酸菌属 |       3 |
| Corynebacterium diphtheriae |       1717 | 白喉棒状杆菌      |      24 |
| Streptomyces hygroscopicus  |       1912 | 吸水链霉菌       |       4 |

## MTBC: assembly

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
        AND ( genus_id in (1763)
            OR species_id in (78331, 1844474, 1305738, 1773, \
                33910, 1717, 1912) )
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND taxonomy_id != species_id               # no strain ID
        AND ( genus_id in (1763)
            OR species_id in (78331, 1844474, 1305738, 1773, \
                33910, 1717, 1912) )
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
    tsv-filter -H --str-ne 3:"Mycobacterium tuberculosis" |
    tsv-filter -H --str-ne 3:"Mycobacterium avium" \
    > raw.1.assembly.tsv

cat raw.0.assembly.tsv |
    tsv-filter -H --str-eq 3:"Mycobacterium tuberculosis" |
    tsv-filter -H --str-in-fld 1:"_variant_" \
    > raw.2.assembly.tsv

# Too many strains
cat raw.0.assembly.tsv |
    tsv-filter -H --str-eq 3:"Mycobacterium tuberculosis" |
    tsv-filter -H --istr-ne 4:"Contig" |
    tsv-filter -H --istr-ne 4:"Scaffold" \
    > raw.3.assembly.tsv

# Poor assembled GA II and IonTorrent reads
cat raw.0.assembly.tsv |
    tsv-filter -H --str-eq 3:"Mycobacterium avium" |
    tsv-filter -H --istr-ne 4:"Contig" |
    tsv-filter -H --istr-ne 4:"Scaffold" \
    > raw.4.assembly.tsv

tsv-append -H raw.{1,2,3,4}.assembly.tsv |
    tsv-uniq -H |
    keep-header -- sort \
    > ${RANK_NAME}.assembly.tsv

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

## Count strains

```bash
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
#  265 ASSEMBLY/n50.pass.csv
#  287 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/MTBC.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/MTBC.assembly.pass.csv

wc -l ASSEMBLY/MTBC.assembly*csv
#   286 ASSEMBLY/MTBC.assembly.collect.csv
#   264 ASSEMBLY/MTBC.assembly.pass.csv

```

```bash
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

#Amycolatopsis   1       3
#Corynebacterium 1       24
#Mycobacterium   19      232
#Streptomyces    1       4

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

```bash
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

```bash
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

```bash
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

```bash
rsync -avP \
    ~/data/alignment/MTBC/ \
    wangq@202.119.37.251:data/alignment/MTBC

# rsync -avP wangq@202.119.37.251:data/alignment/MTBC/ ~/data/alignment/MTBC

```

* `--perseq` for Chromosome-level assemblies and targets

```bash
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

```bash
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

```bash
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

