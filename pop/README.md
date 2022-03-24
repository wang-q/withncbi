# Build alignments across a eukaryotic taxonomy rank

Genus *Trichoderma* as an example.

[TOC levels=1-3]: # ""

- [Build alignments across a eukaryotic taxonomy rank](#build-alignments-across-a-eukaryotic-taxonomy-rank)
    * [Preparations](#preparations)
    * [Strain info](#strain-info)
        + [List all ranks](#list-all-ranks)
        + [Species with assemblies](#species-with-assemblies)
    * [Trichoderma: assembly](#trichoderma-assembly)
    * [Filter strains by N50](#filter-strains-by-n50)
    * [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
    * [Groups and targets](#groups-and-targets)
    * [Prepare sequences for `egaz`](#prepare-sequences-for-egaz)
    * [Generate alignments](#generate-alignments)

## Preparations

* Install `nwr` and create a local taxonomy and assembly database.

```shell
brew install wang-q/tap/nwr
brew install sqlite

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

```

## Strain info

* [Trichoderma](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=5543)
* [Entrez records](http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5543)

### List all ranks

There are no noteworthy classification ranks other than species.

```shell
mkdir -p ~/data/alignment/Trichoderma
cd ~/data/alignment/Trichoderma

nwr member Trichoderma |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr lineage Trichoderma |
    tsv-filter --str-ne 1:clade |
    sed -n '/kingdom\tFungi/,$p' |
    (echo -e '#rank\tsci_name\ttax_id' && cat) |
    mlr --itsv --omd cat

```

| rank     | count |
|----------|------:|
| genus    |     1 |
| species  |   420 |
| no rank  |     1 |
| varietas |     2 |
| strain   |    14 |
| forma    |     2 |

| #rank      | sci_name          | tax_id |
|------------|-------------------|--------|
| kingdom    | Fungi             | 4751   |
| subkingdom | Dikarya           | 451864 |
| phylum     | Ascomycota        | 4890   |
| subphylum  | Pezizomycotina    | 147538 |
| class      | Sordariomycetes   | 147550 |
| subclass   | Hypocreomycetidae | 222543 |
| order      | Hypocreales       | 5125   |
| family     | Hypocreaceae      | 5129   |
| genus      | Trichoderma       | 5543   |

### Species with assemblies

Check also the family Hypocreaceae for outgroups.

```shell
cd ~/data/alignment/Trichoderma

SPECIES=$(
    nwr member Hypocreaceae -r species |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        sed '1d' |
        cut -f 1 |
        sort
)

for S in $SPECIES; do
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

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_genbank.sqlite
    )

    if [[ ${GB} -gt 0 ]]; then
        echo -e "$S\t$GB\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e '#tax_id\tspecies\tGB\tCHR' && cat) \
    > species.count.tsv

cat species.count.tsv |
    tsv-filter -H --ge GB:1 |
    sed 's/Trichoderma /T. /g' |
    mlr --itsv --omd cat

```

| #tax_id | species                | GB  | CHR |
|---------|------------------------|-----|-----|
| 51453   | T. reesei              | 11  | 7   |
| 101201  | T. asperellum          | 11  | 2   |
| 5544    | T. harzianum           | 9   | 1   |
| 63577   | T. atroviride          | 7   | 1   |
| 29875   | T. virens              | 6   | 2   |
| 150374  | Escovopsis weberi      | 2   | 0   |
| 1567482 | T. afroharzianum       | 2   | 0   |
| 398673  | T. gamsii              | 2   | 0   |
| 49224   | T. hamatum             | 2   | 0   |
| 5548    | T. longibrachiatum     | 2   | 0   |
| 654480  | T. cornu-damae         | 1   | 1   |
| 1491008 | T. semiorbis           | 1   | 1   |
| 1491479 | T. simmonsii           | 1   | 1   |
| 767780  | Cladobotryum protrusum | 1   | 0   |
| 2060699 | Hypomyces perniciosus  | 1   | 0   |
| 5132    | Hypomyces rosellus     | 1   | 0   |
| 490622  | T. arundinaceum        | 1   | 0   |
| 702382  | T. asperelloides       | 1   | 0   |
| 1491457 | T. atrobrunneum        | 1   | 0   |
| 247546  | T. brevicompactum      | 1   | 0   |
| 2034171 | T. brevicrassum        | 1   | 0   |
| 58853   | T. citrinoviride       | 1   | 0   |
| 202914  | T. erinaceum           | 1   | 0   |
| 1195189 | T. gracile             | 1   | 0   |
| 1491466 | T. guizhouense         | 1   | 0   |
| 97093   | T. koningii            | 1   | 0   |
| 337941  | T. koningiopsis        | 1   | 0   |
| 1567552 | T. lentiforme          | 1   | 0   |
| 1491472 | T. lixii               | 1   | 0   |
| 1497375 | T. oligosporum         | 1   | 0   |
| 858221  | T. parareesei          | 1   | 0   |
| 500994  | T. pleuroti            | 1   | 0   |
| 5547    | T. viride              | 1   | 0   |

## Download all assemblies

[WGS](https://www.ncbi.nlm.nih.gov/Traces/wgs/?view=wgs&search=Trichoderma) is now useless and can
be ignored.

```shell
cd ~/data/alignment/Trichoderma

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND (
            (genus IN ('Trichoderma'))
            OR
            (species IN ('Escovopsis weberi', 'Cladobotryum protrusum', 'Hypomyces perniciosus', 'Hypomyces rosellus'))
        )
        AND species NOT LIKE '% sp.%'
    " |
    sqlite3 -tabs ~/.nwr/ar_genbank.sqlite \
    > raw.tsv

cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
    perl ~/Scripts/withncbi/taxon/abbr_name.pl -c "1,2,3" -s '\t' -m 3 --shortsub |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"," -e '
        BEGIN{my %seen};
        /^#/ and print and next;
        /^organism_name/i and next;
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
        $seen{$F[5]} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $F[5], $F[3], $F[1], $F[4];
        ' |
    keep-header -- sort -k3,3 -k1,1 \
    > Trichoderma.assembly.tsv

# find potential duplicated strains or assemblies
cat Trichoderma.assembly.tsv |
    tsv-uniq -f 1 --repeated

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Trichoderma.assembly.tsv
# cp Trichoderma.assembly.tsv ~/Scripts/withncbi/pop

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

Information of assemblies are collected from *_assembly_report.txt *after* downloading.

**Note**: `*_assembly_report.txt` have `CRLF` at the end of the line.

```shell
cd ~/data/alignment/Trichoderma

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Trichoderma.assembly.tsv \
    -o ASSEMBLY

# Remove dirs not in the list
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    tr "/" "\t" |
    cut -f 2 |
    tsv-join --exclude -k 1 -f ASSEMBLY/rsync.tsv -d 1 |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo Remove {}
        rm -fr ASSEMBLY/{}
    '

# Run
proxychains4 bash ASSEMBLY/Pseudomonas.assembly.rsync.sh

bash ASSEMBLY/Pseudomonas.assembly.collect.sh

# md5
cat ASSEMBLY/rsync.tsv |
    tsv-select -f 1 |
    parallel -j 4 --keep-order '
        echo "==> {}"
        cd ASSEMBLY/{}
        md5sum --check md5checksums.txt
    ' |
    grep -v ": OK"

```

## Count and group strains

```shell
cd ~/data/alignment/Trichoderma

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
#  66 ASSEMBLY/n50.pass.csv
#  78 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/Trichoderma.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv \
    > ASSEMBLY/Trichoderma.assembly.pass.csv

wc -l ASSEMBLY/Trichoderma.assembly*csv
#   78 ASSEMBLY/Trichoderma.assembly.collect.csv
#   66 ASSEMBLY/Trichoderma.assembly.pass.csv

```


* strains

```shell
cd ~/data/alignment/Trichoderma

# list strains
mkdir -p taxon

rm taxon/*
cat ASSEMBLY/Trichoderma.assembly.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,2,3 |
    nwr append stdin -c 3 -r species -r genus -r family -r order |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 1 '
        if [[ "{#}" -eq "1" ]]; then
            rm strains.lst
            rm genus.tmp
            rm species.tmp
        fi

        echo {1} >> strains.lst

        echo {5} >> genus.tmp
        echo {1} >> taxon/{5}

        echo {4} >> species.tmp

        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\n" {1} {2} {3} {4} {5} {6} {7}
    ' \
    > strains.taxon.tsv

cat genus.tmp | tsv-uniq > genus.lst
cat species.tmp | tsv-uniq > species.lst

# Omit strains without protein annotations
for STRAIN in $(cat strains.lst); do
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
        echo ${STRAIN}
    fi
    if ! compgen -G "ASSEMBLY/${STRAIN}/*_cds_from_genomic.fna.gz" > /dev/null; then
        echo ${STRAIN}
    fi
done |
    tsv-uniq \
    > omit.lst

rm *.tmp

```

## NCBI taxonomy

Done by `bp_taxonomy2tree.pl` from BioPerl.

```shell
mkdir -p ~/data/alignment/Trichoderma/tree
cd ~/data/alignment/Trichoderma/tree

bp_taxonomy2tree.pl -e \
    $(
        cat ../species.lst |
            tr " " "_" |
            parallel echo '-s {}'
    ) |
    sed 's/Trichoderma/T/g' \
    > ncbi.nwk

nw_display -s -b 'visibility:hidden' -w 600 -v 30 ncbi.nwk |
    rsvg-convert -o Trichoderma.ncbi.png

```

## Raw phylogenetic tree by MinHash

```shell
mkdir -p ~/data/alignment/Trichoderma/mash
cd ~/data/alignment/Trichoderma/mash

for strain in $(cat ../strains.lst ); do
    2>&1 echo "==> ${strain}"

    if [[ -e ${strain}.msh ]]; then
        continue
    fi

    find ../ASSEMBLY/${strain} -name "*_genomic.fna.gz" |
        grep -v "_from_" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 8 - -I "${strain}" -o ${strain}
done

mash triangle -E -p 8 -l <(
    cat ../strains.lst | parallel echo "{}.msh"
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

```

### Tweak the mash tree

```shell
mkdir -p ~/data/alignment/Trichoderma/tree
cd ~/data/alignment/Trichoderma/tree

nw_reroot ../mash/tree.nwk C_pro_GCA_004303015_1 H_per_GCA_008477525_1 |
    nw_order -c n - \
    > mash.reroot.newick

# rank::col
ARRAY=(
#    'order::7'
#    'family::6'
    'genus::5'
    'species::4'
)

rm mash.condensed.map
CUR_TREE=mash.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick mash.${GROUP_NAME}.newick
    cat condense.map >> mash.condensed.map

    CUR_TREE=mash.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 mash.species.newick |
    rsvg-convert -o Trichoderma.mash.png

# cp Trichoderma.mash.png ~/Scripts/withncbi/image/Trichoderma.png

```

![Trichoderma.png](../image/Trichoderma.png)

## Groups and targets

Review `ASSEMBLY/Trichoderma.assembly.pass.csv` and `mash/groups.tsv`.

Create `ARRAY` manually with a format `group::target`.

Target criteria:

* Prefer Sander sequenced assemblies
* RefSeq_category with `Representative Genome`
* Assembly_level with `Complete Genome` or `Chromosome`

```shell
mkdir -p ~/data/alignment/Trichoderma/taxon
cd ~/data/alignment/Trichoderma/taxon

cp ../mash/tree.nwk .
cp ../mash/groups.tsv .

echo -e "#Serial\tGroup\tCount\tTarget" > group_target.tsv

# groups accroding `groups.tsv`
ARRAY=(
    'C_E_H::H_ros_GCA_011799845_1'
    'T_har_vire::T_vire_Gv29_8_GCA_000170995_2'
    'T_asperellum_atrov::T_atrov_IMI_206040_GCA_000171015_2'
    'T_lon_ree::T_ree_QM6a_GCA_000167675_2'
)

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
    'Trichoderma::T_ree_QM6a_GCA_000167675_2'
    'Trichoderma_reesei::T_ree_QM6a_GCA_000167675_2'
    'Trichoderma_asperellum::T_asperellum_CBS_433_97_GCA_003025105_1'
    'Trichoderma_harzianum::T_har_CBS_226_95_GCA_003025095_1'
    'Trichoderma_atroviride::T_atrov_IMI_206040_GCA_000171015_2'
    'Trichoderma_virens::T_vire_Gv29_8_GCA_000170995_2'
)

SERIAL=100
for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    TARGET_NAME="${item##*::}"

    SERIAL=$((SERIAL + 1))
    GROUP_NAME_2=$(echo $GROUP_NAME | tr "_" " ")

    if [ "$GROUP_NAME" = "Trichoderma" ]; then
        cat ../ASSEMBLY/Trichoderma.assembly.pass.csv |
            tsv-filter -H -d"," --not-blank RefSeq_category |
            sed '1d' |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
        echo "C_pro_GCA_004303015_1" >> ${GROUP_NAME}
        echo "E_web_GCA_003055145_1" >> ${GROUP_NAME}
        echo "H_per_GCA_008477525_1" >> ${GROUP_NAME}
        echo "H_ros_GCA_011799845_1" >> ${GROUP_NAME}
    else
        cat ../ASSEMBLY/Trichoderma.assembly.pass.csv |
            cut -d"," -f 1,2 |
            grep "${GROUP_NAME_2}" |
            cut -d"," -f 1 \
            > ${GROUP_NAME}
    fi

    COUNT=$(cat ${GROUP_NAME} | wc -l )

    echo -e "${SERIAL}\t${GROUP_NAME}\t${COUNT}\t${TARGET_NAME}" >> group_target.tsv

done

mlr --itsv --omd cat group_target.tsv

cat ../ASSEMBLY/Trichoderma.assembly.pass.csv |
    tsv-filter -H -d, --str-eq Assembly_level:"Complete Genome" |
    tsv-select -H -d, -f name \
    > complete-genome.lst


```

| #Serial | Group                  | Count | Target                                  |
|---------|------------------------|-------|-----------------------------------------|
| 1       | C_E_H                  | 6     | H_ros_GCA_011799845_1                   |
| 2       | T_har_vire             | 20    | T_vire_Gv29_8_GCA_000170995_2           |
| 3       | T_asperellum_atrov     | 23    | T_atrov_IMI_206040_GCA_000171015_2      |
| 4       | T_lon_ree              | 16    | T_ree_QM6a_GCA_000167675_2              |
| 101     | Trichoderma            | 31    | T_ree_QM6a_GCA_000167675_2              |
| 102     | Trichoderma_reesei     | 10    | T_ree_QM6a_GCA_000167675_2              |
| 103     | Trichoderma_asperellum | 10    | T_asperellum_CBS_433_97_GCA_003025105_1 |
| 104     | Trichoderma_harzianum  | 7     | T_har_CBS_226_95_GCA_003025095_1        |
| 105     | Trichoderma_atroviride | 7     | T_atrov_IMI_206040_GCA_000171015_2      |
| 106     | Trichoderma_virens     | 6     | T_vire_Gv29_8_GCA_000170995_2           |

## Prepare sequences for `egaz`

* `--perseq` for Chromosome-level assemblies and targets
    * means split fasta by names, target or good assembles should set it
* `--species Fungi` specify the species or clade of this group for RepeatMasker

```shell
cd ~/data/alignment/Trichoderma/

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    $( cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 | parallel -j 1 echo " --perseq {} " ) \
    $( cat taxon/complete-genome.lst | parallel -j 1 echo " --perseq {} " ) \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Fungi --parallel 24"

bash GENOMES/0_prep.sh

# gff
for n in $(cat taxon/group_target.tsv | sed -e '1d' | cut -f 4 ) \
    $( cat taxon/complete-genome.lst ) \
    ; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"

    gzip -dc ${FILE_GFF} > GENOMES/${n}/chr.gff
done

```

## Generate alignments

```shell
cd ~/data/alignment/Trichoderma/

cat taxon/group_target.tsv |
    sed -e '1d' |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo -e "==> Group: [{2}]\tTarget: [{4}]\n"

        egaz template \
            GENOMES/{4} \
            $(cat taxon/{2} | grep -v -x "{4}" | xargs -I[] echo "GENOMES/[]") \
            --multi -o groups/{2}/ \
            --tree taxon/tree.nwk \
            --parallel 8 -v

        bash groups/{2}/1_pair.sh
        bash groups/{2}/3_multi.sh
    '

# clean
find groups -mindepth 1 -maxdepth 3 -type d -name "*_raw" | parallel -r rm -fr
find groups -mindepth 1 -maxdepth 3 -type d -name "*_fasta" | parallel -r rm -fr
find . -mindepth 1 -maxdepth 3 -type f -name "output.*" | parallel -r rm

```

Use Tatr_IMI_206040 as target

* No results for Tatr_IMI_206040vsTkon_JCM_1883

Tatr_IMI_206040;qs=Tatr_XS2015,,

```bash
cd ~/data/alignment/trichoderma

egaz template \
    GENOMES/Tatr_IMI_206040 \
    GENOMES/Tree_QM6a \
    GENOMES/Tvir_Gv29_8 \
    --multi -o multi/ \
    --multiname sanger --order \
    --parallel 8 -v

bash multi/1_pair.sh
bash multi/3_multi.sh

egaz template \
    GENOMES/Tatr_IMI_206040 \
    $(find GENOMES -maxdepth 1 -mindepth 1 -type d -path "*/????*" | grep -v "Tatr_IMI_206040"| grep -v "Tkon_JCM_1883") \
    --multi -o multi/ \
    --tree mash/tree.nwk \
    --parallel 8 -v

bash multi/1_pair.sh
bash multi/3_multi.sh

```
