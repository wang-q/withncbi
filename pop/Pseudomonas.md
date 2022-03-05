# *Pseudomonas* HGT

- [*Pseudomonas* HGT](#pseudomonas-hgt)
    * [Software](#software)
    * [Strain info](#strain-info)
        + [List all ranks](#list-all-ranks)
        + [Species with assemblies](#species-with-assemblies)
        + [Outgroups](#outgroups)
    * [Download all assemblies](#download-all-assemblies)
    * [Count and group strains](#count-and-group-strains)
    * [NCBI taxonomy](#ncbi-taxonomy)
    * [Raw phylogenetic tree by MinHash](#raw-phylogenetic-tree-by-minhash)
        + [Tweak the mash tree](#tweak-the-mash-tree)
    * [Collect proteins](#collect-proteins)
        + [`all.pro.fa`](#allprofa)
        + [`all.replace.fa`](#allreplacefa)
        + [`all.info.tsv`](#allinfotsv)
    * [Phylogenetics with 40 single-copy genes and RNase_R](#phylogenetics-with-40-single-copy-genes-and-rnase_r)
        + [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
        + [Create a valid marker gene list](#create-a-valid-marker-gene-list)
        + [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
        + [Tweak the concat tree](#tweak-the-concat-tree)

The genus Pseudomonas includes the conditionally pathogenic bacteria Pseudomonas aeruginosa, plant
pathogens, plant beneficial bacteria, and soil bacteria. Microorganisms of Pseudomonas are extremely
rich in metabolic diversity, and it is thought that this diversity also allows them to survive in a
very wide range of ecological niches.

The metabolism of carbohydrates is a fundamental biochemical process that ensures a continuous
supply of energy to living cells.

The biodegradation of RNA is also an important part of metabolism, which is accomplished by the
degradosomes. However, the diversity of degradosomes in different environments has not been fully
investigated.

According to a recent [paper](https://doi.org/10.1128/mSystems.00543-20), there are some order-level
changes in Gammaproteobacteria. We include both old and new orders.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

## Software

* Install `nwr` and create a local taxonomy and assembly database.

```shell
brew install wang-q/tap/nwr # 0.5.5 or above
brew install sqlite         # 3.34 or above

nwr download
nwr txdb

nwr ardb
nwr ardb --genbank

```

* Other packages

```shell
brew install hmmer
brew install samtools
brew install librsvg

brew install brewsci/bio/muscle
brew install brewsci/bio/fasttree
brew install brewsci/bio/easel
brew install brewsci/bio/newick-utils
brew install brewsci/bio/trimal

brew install datamash
brew install miller
brew install wang-q/tap/tsv-utils

brew install jq
brew install pup

```

## Strain info

* [Pseudomonas](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=286)
* [Acinetobacter](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=469)

### List all ranks

```shell
mkdir -p ~/data/Pseudomonas
cd ~/data/Pseudomonas

nwr member Pseudomonas |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Acinetobacter |
    grep -v " sp." |
    tsv-summarize -H -g 3 --count |
    mlr --itsv --omd cat

nwr member Pseudomonas Acinetobacter -r "species group" -r "species subgroup" |
    tsv-select -f 1-3 |
    keep-header -- tsv-sort -k3,3 -k2,2 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species          | 405   |
| strain           | 747   |
| subspecies       | 12    |
| no rank          | 120   |
| species group    | 6     |
| species subgroup | 5     |
| isolate          | 1     |

| rank             | count |
|------------------|-------|
| genus            | 1     |
| species group    | 2     |
| species subgroup | 3     |
| species          | 113   |
| strain           | 1110  |
| no rank          | 2     |
| subspecies       | 1     |
| isolate          | 2     |

| #tax_id | sci_name                                 | rank             |
|---------|------------------------------------------|------------------|
| 909768  | A. calcoaceticus/baumannii complex       | species group    |
| 2839056 | A. Taxon 24                              | species group    |
| 136841  | P. aeruginosa group                      | species group    |
| 136842  | P. chlororaphis group                    | species group    |
| 136843  | P. fluorescens group                     | species group    |
| 136845  | P. putida group                          | species group    |
| 136846  | P. stutzeri group                        | species group    |
| 136849  | P. syringae group                        | species group    |
| 2839060 | A. Taxon 24C                             | species subgroup |
| 2839057 | A. Taxon 24D                             | species subgroup |
| 2839061 | A. Taxon 24E                             | species subgroup |
| 627141  | P. nitroreducens/multiresinivorans group | species subgroup |
| 1232139 | P. oleovorans/pseudoalcaligenes group    | species subgroup |
| 578833  | P. stutzeri subgroup                     | species subgroup |
| 251695  | P. syringae group genomosp. 1            | species subgroup |
| 251698  | P. syringae group genomosp. 2            | species subgroup |

### Species with assemblies

Check also the order Pseudomonadales.

* Old ones: Cellvibrionales, Oceanospirillales, Pseudomonadales, and Alteromonadales
* New ones: Moraxellales, Kangiellales, and Pseudomonadales

```shell
cd ~/data/Pseudomonas

SPECIES=$(
    nwr member -r species \
        Cellvibrionales Oceanospirillales Alteromonadales \
        Moraxellales Kangiellales Pseudomonadales |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        grep -v " sp." |
        grep -v -E "\bbacterium\b" |
        grep -v -E "\bsymbiont\b" |
        sed '1d' |
        cut -f 1 |
        sort |
        uniq
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

    CHR=$(
        echo "
            SELECT
                COUNT(*)
            FROM ar
            WHERE 1=1
                AND species_id = $S
                AND assembly_level IN ('Complete Genome', 'Chromosome')
            " |
            sqlite3 -tabs ~/.nwr/ar_refseq.sqlite
    )

    if [[ ${RS} -gt 0 ]]; then
        echo -e "$S\t$RS\t$CHR"
    fi
done |
    nwr append stdin |
    tsv-select -f 1,4,2-3 |
    tsv-sort -k3,3nr -k4,4nr -k2,2 |
    (echo -e '#tax_id\tspecies\tRS\tCHR' && cat) \
    > species.count.tsv

cat species.count.tsv |
    tsv-filter -H --ge CHR:5 |
    tsv-filter -H --invert --str-in-fld species:Pseudomonas --lt RS:30 |
    tsv-filter -H --invert --str-in-fld species:Acinetobacter --lt RS:30 |
    sed 's/Pseudomonas /P. /g' |
    sed 's/Acinetobacter /A. /g' |
    mlr --itsv --omd cat

```

| #tax_id | species                         | RS   | CHR |
|---------|---------------------------------|------|-----|
| 287     | P. aeruginosa                   | 6310 | 475 |
| 470     | A. baumannii                    | 5933 | 332 |
| 33069   | P. viridiflava                  | 1536 | 7   |
| 317     | P. syringae                     | 540  | 41  |
| 48296   | A. pittii                       | 315  | 32  |
| 294     | P. fluorescens                  | 257  | 39  |
| 480     | Moraxella catarrhalis           | 210  | 16  |
| 303     | P. putida                       | 189  | 49  |
| 106654  | A. nosocomialis                 | 159  | 11  |
| 316     | P. stutzeri                     | 133  | 30  |
| 29438   | P. savastanoi                   | 116  | 5   |
| 38313   | Shewanella algae                | 108  | 22  |
| 587753  | P. chlororaphis                 | 99   | 60  |
| 756892  | A. indicus                      | 91   | 19  |
| 47877   | P. amygdali                     | 86   | 8   |
| 380021  | P. protegens                    | 73   | 23  |
| 40215   | A. junii                        | 65   | 9   |
| 1530123 | A. seifertii                    | 58   | 25  |
| 29430   | A. haemolyticus                 | 55   | 14  |
| 76759   | P. monteilii                    | 46   | 9   |
| 40214   | A. johnsonii                    | 43   | 19  |
| 296     | P. fragi                        | 38   | 6   |
| 28090   | A. lwoffii                      | 31   | 11  |
| 43657   | Pseudoalteromonas luteoviolacea | 25   | 5   |
| 28108   | Alteromonas macleodii           | 24   | 9   |
| 34062   | Moraxella osloensis             | 21   | 10  |
| 43662   | Pseudoalteromonas piscicida     | 18   | 6   |
| 314275  | Alteromonas mediterranea        | 16   | 16  |
| 24      | Shewanella putrefaciens         | 15   | 9   |
| 62322   | Shewanella baltica              | 14   | 11  |
| 386891  | Moraxella bovoculi              | 9    | 7   |
| 1697053 | Thiopseudomonas alkaliphila     | 7    | 7   |

### Outgroups

Use these model organisms as outgroups.

```shell
cd ~/data/Pseudomonas

GENUS=$(
    nwr member Bacteria -r genus |
        grep -v -i "Candidatus " |
        grep -v -i "candidate " |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
.headers ON

    SELECT
        *
    FROM ar
    WHERE 1=1
        AND genus_id IN ($GENUS)
        AND refseq_category IN ('reference genome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite \
    > reference.tsv

cat reference.tsv |
    sed '1s/^/#/' |
    nwr append stdin -r phylum -r class |
    tsv-select -H -f 1,2,phylum,class |
    parallel --col-sep "\t" -j 1 '
        if [[ "{3}" != "Proteobacteria" ]]; then
            printf "%s\t%s\t%s\n" {1} {2} {3}
        else
            printf "%s\t%s\t%s\n" {1} {2} {4}
        fi
    ' |
    mlr --itsv --omd cat

```

| #tax_id | organism_name                                                    | phylum                |
|---------|------------------------------------------------------------------|-----------------------|
| 565050  | Caulobacter vibrioides NA1000                                    | Alphaproteobacteria   |
| 192222  | Campylobacter jejuni subsp. jejuni NCTC 11168 = ATCC 700819      | Epsilonproteobacteria |
| 208964  | Pseudomonas aeruginosa PAO1                                      | Gammaproteobacteria   |
| 871585  | Acinetobacter pittii PHEA-2                                      | Gammaproteobacteria   |
| 511145  | Escherichia coli str. K-12 substr. MG1655                        | Gammaproteobacteria   |
| 386585  | Escherichia coli O157:H7 str. Sakai                              | Gammaproteobacteria   |
| 1125630 | Klebsiella pneumoniae subsp. pneumoniae HS11286                  | Gammaproteobacteria   |
| 99287   | Salmonella enterica subsp. enterica serovar Typhimurium str. LT2 | Gammaproteobacteria   |
| 198214  | Shigella flexneri 2a str. 301                                    | Gammaproteobacteria   |
| 227377  | Coxiella burnetii RSA 493                                        | Gammaproteobacteria   |
| 272561  | Chlamydia trachomatis D/UW-3/CX                                  | Chlamydiae            |
| 93061   | Staphylococcus aureus subsp. aureus NCTC 8325                    | Firmicutes            |
| 224308  | Bacillus subtilis subsp. subtilis str. 168                       | Firmicutes            |
| 169963  | Listeria monocytogenes EGD-e                                     | Firmicutes            |
| 83332   | Mycobacterium tuberculosis H37Rv                                 | Actinobacteria        |

## Download all assemblies

Species with 2 or more genomes were retained.

```shell
cd ~/data/Pseudomonas

SPECIES=$(
    cat species.count.tsv |
        tsv-filter -H --ge CHR:2 |
        sed '1d' |
        cut -f 1 |
        tr "\n" "," |
        sed 's/,$//'
)

echo "
    SELECT
        organism_name || ' ' || assembly_accession AS name,
        species, genus, ftp_path, assembly_level
    FROM ar
    WHERE 1=1
        AND species_id IN ($SPECIES)
        AND species NOT LIKE '% sp.%'
        AND organism_name NOT LIKE '% sp.%'
        AND assembly_level IN ('Complete Genome', 'Chromosome')
    " |
    sqlite3 -tabs ~/.nwr/ar_refseq.sqlite |
    tsv-filter --invert --str-eq 2:"Pseudomonas aeruginosa" --str-eq 5:"Chromosome" |
    tsv-filter --invert --str-eq 2:"Acinetobacter baumannii" --str-eq 5:"Chromosome" \
    > raw.tsv

cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level \
    >> raw.tsv

cat raw.tsv |
    grep -v '^#' |
    tsv-uniq |
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
    > Pseudomonas.assembly.tsv

# find potential duplicated strains or assemblies
cat Pseudomonas.assembly.tsv |
    tsv-uniq -f 1 --repeated

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim Pseudomonas.assembly.tsv
# cp Pseudomonas.assembly.tsv ~/Scripts/withncbi/pop

# Comment out unneeded strains

# Cleaning
rm raw*.*sv

```

```shell
cd ~/data/Pseudomonas

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Pseudomonas.assembly.tsv \
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

* Check N50 of assemblies

* Pseudom_flu_GCF_900636635_1 is a strain of P. aeruginosa

```shell
cd ~/data/Pseudomonas

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
#  1522 ASSEMBLY/n50.pass.csv
#  1522 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/Pseudomonas.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_900636635 \
    > ASSEMBLY/Pseudomonas.assembly.pass.csv

wc -l ASSEMBLY/Pseudomonas.assembly*csv
#  1522 ASSEMBLY/Pseudomonas.assembly.collect.csv
#  1521 ASSEMBLY/Pseudomonas.assembly.pass.csv

```

* Genus

```shell
cd ~/data/Pseudomonas

# Group by genus
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tsv-select -d, -f 3 |
    tsv-uniq |
    nwr append stdin -r genus |
    tsv-select -f 2 |
    tsv-uniq \
    > genus.lst

cat genus.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        n_species=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus -r species |
            grep {} |
            tsv-select -f 1,3 |
            tsv-uniq |
            wc -l)

        n_strains=$(cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
            sed "1d" |
            tsv-select -d, -f 3 |
            nwr append stdin -r genus |
            grep {} |
            wc -l)

        printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' |
    nwr append stdin --id |
    tsv-select -f 5,4,2,3 |
    tsv-sort -k2,2 |
    (echo -e '#tax_id\tgenus\t#species\t#strains' && cat) |
    mlr --itsv --omd cat

```

| #tax_id | genus             | #species | #strains |
|---------|-------------------|----------|----------|
| 469     | Acinetobacter     | 43       | 484      |
| 226     | Alteromonas       | 16       | 31       |
| 352     | Azotobacter       | 5        | 6        |
| 1386    | Bacillus          | 1        | 1        |
| 194     | Campylobacter     | 1        | 1        |
| 75      | Caulobacter       | 1        | 1        |
| 10      | Cellvibrio        | 2        | 4        |
| 810     | Chlamydia         | 1        | 1        |
| 776     | Coxiella          | 1        | 1        |
| 561     | Escherichia       | 2        | 2        |
| 2745    | Halomonas         | 5        | 13       |
| 135575  | Idiomarina        | 2        | 2        |
| 570     | Klebsiella        | 1        | 1        |
| 1637    | Listeria          | 1        | 1        |
| 2742    | Marinobacter      | 5        | 6        |
| 28253   | Marinomonas       | 1        | 2        |
| 475     | Moraxella         | 5        | 35       |
| 1763    | Mycobacterium     | 1        | 1        |
| 53246   | Pseudoalteromonas | 18       | 33       |
| 286     | Pseudomonas       | 172      | 829      |
| 497     | Psychrobacter     | 2        | 2        |
| 590     | Salmonella        | 1        | 1        |
| 22      | Shewanella        | 15       | 50       |
| 620     | Shigella          | 1        | 1        |
| 1279    | Staphylococcus    | 1        | 1        |
| 187492  | Thalassolituus    | 3        | 3        |
| 1654787 | Thiopseudomonas   | 1        | 7        |

* strains

```shell
cd ~/data/Pseudomonas

# list strains
mkdir -p taxon

rm taxon/*
cat ASSEMBLY/Pseudomonas.assembly.pass.csv |
    sed -e '1d' |
    tr "," "\t" |
    tsv-select -f 1,2,3 |
    nwr append stdin -c 3 -r genus -r "species group" -r species |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 1 '
        if [[ "{#}" -eq "1" ]]; then
            rm strains.lst
            rm genus.tmp
            rm species_group.tmp
            rm species.tmp
        fi

        echo {1} >> strains.lst

        echo {4} >> genus.tmp
        echo {1} >> taxon/{4}

        echo {5} >> species_group.tmp
        echo {6} >> species.tmp

        printf "%s\t%s\t%d\t%s\t%s\t%s\n" {1} {2} {3} {4} {5} {6}
    ' \
    > strains.taxon.tsv

cat genus.tmp | tsv-uniq > genus.lst
cat species_group.tmp | tsv-uniq > species_group.lst
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
# All OK

rm *.tmp

```

## NCBI taxonomy

Done by `bp_taxonomy2tree.pl` from BioPerl.

```shell
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

bp_taxonomy2tree.pl -e \
    $(
        cat ../species.lst |
            tr " " "_" |
            parallel echo '-s {}'
    ) |
    sed 's/Pseudomonas/P/g' \
    > ncbi.nwk

nw_display -s -b 'visibility:hidden' -w 600 -v 30 ncbi.nwk |
    rsvg-convert -o Pseudomonas.ncbi.png

```

## Raw phylogenetic tree by MinHash

```shell
mkdir -p ~/data/Pseudomonas/mash
cd ~/data/Pseudomonas/mash

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
mkdir -p ~/data/Pseudomonas/tree
cd ~/data/Pseudomonas/tree

nw_reroot ../mash/tree.nwk B_sub_subtilis_168 St_aur_aureus_NCTC_8325 \
    > mash.reroot.newick

rm mash.condensed.map

# genus
bash ~/Scripts/withncbi/taxon/condense_tree.sh mash.reroot.newick ../strains.taxon.tsv 1 4

mv mash.condense.newick mash.genus.newick
cat mash.condense.map >> mash.condensed.map

# species_group
bash ~/Scripts/withncbi/taxon/condense_tree.sh mash.genus.newick ../strains.taxon.tsv 1 5

mv mash.genus.condense.newick mash.species_group.newick
cat mash.genus.condense.map >> mash.condensed.map

# species
bash ~/Scripts/withncbi/taxon/condense_tree.sh mash.species_group.newick ../strains.taxon.tsv 1 6

mv mash.species_group.condense.newick mash.species.newick
cat mash.species_group.condense.map >> mash.condensed.map

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 mash.species.newick |
    rsvg-convert -o Pseudomonas.mash.png

```

## Collect proteins

### `all.pro.fa`

```shell script
cd ~/data/Pseudomonas

mkdir -p PROTEINS

# 1521
find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l

# 1521
find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l

# 1520
cat strains.lst |
    wc -l

for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz
done \
    > PROTEINS/all.pro.fa

cat PROTEINS/all.pro.fa |
    perl -nl -e '
        BEGIN { our %seen; our $h; }

        if (/^>/) {
            $h = (split(" ", $_))[0];
            $seen{$h}++;
            $_ = $h;
        }
        print if $seen{$h} == 1;
    ' \
    > PROTEINS/all.uniq.fa

# counting proteins
cat PROTEINS/all.pro.fa |
    grep "^>" |
    wc -l
#7182872

cat PROTEINS/all.pro.fa |
    grep "^>" |
    tsv-uniq |
    wc -l
#2487274

# annotations may be different
cat PROTEINS/all.uniq.fa |
    grep "^>" |
    wc -l
#2438570

# ribonuclease
cat PROTEINS/all.pro.fa |
    grep "ribonuclease" |
    grep -v "deoxyribonuclease" |
    perl -nl -e 's/^>\w+\.\d+\s+//g; print' |
    perl -nl -e 's/\s+\[.+?\]$//g; print' |
    perl -nl -e 's/MULTISPECIES: //g; print' |
    sort |
    uniq -c |
    sort -nr

```

### `all.replace.fa`

```shell
cd ~/data/Pseudomonas

rm PROTEINS/all.strain.tsv
for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        STRAIN=${STRAIN} perl -nl -e '
            $n = $_;
            $s = $n;
            $s =~ s/\.\d+//;
            printf qq{%s\t%s_%s\t%s\n}, $n, $ENV{STRAIN}, $s, $ENV{STRAIN};
        ' \
    > PROTEINS/${STRAIN}.replace.tsv

    cut -f 2,3 PROTEINS/${STRAIN}.replace.tsv >> PROTEINS/all.strain.tsv

    faops replace -s ASSEMBLY/${STRAIN}/*_protein.faa.gz <(cut -f 1,2 PROTEINS/${STRAIN}.replace.tsv) stdout

    rm PROTEINS/${STRAIN}.replace.tsv
done \
    > PROTEINS/all.replace.fa

cat PROTEINS/all.replace.fa |
    grep "^>" |
    wc -l
#7182872

(echo -e "#name\tstrain" && cat PROTEINS/all.strain.tsv)  \
    > temp &&
    mv temp PROTEINS/all.strain.tsv

faops size PROTEINS/all.replace.fa > PROTEINS/all.replace.sizes

(echo -e "#name\tsize" && cat PROTEINS/all.replace.sizes) > PROTEINS/all.size.tsv

rm PROTEINS/all.replace.sizes

```

### `all.info.tsv`

```shell
cd ~/data/Pseudomonas

for STRAIN in $(cat strains.lst); do
    gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
        grep "^>" |
        sed "s/^>//" |
        perl -nl -e '/\[.+\[/ and s/\[/\(/; print' |
        perl -nl -e '/\].+\]/ and s/\]/\)/; print' |
        perl -nl -e 's/\s+\[.+?\]$//g; print' |
        perl -nl -e 's/MULTISPECIES: //g; print' |
        STRAIN=${STRAIN} perl -nl -e '
            /^(\w+)\.\d+\s+(.+)$/ or next;
            printf qq{%s_%s\t%s\n}, $ENV{STRAIN}, $1, $2;
        '
done \
    > PROTEINS/all.annotation.tsv

cat PROTEINS/all.annotation.tsv |
    wc -l
#7182872

(echo -e "#name\tannotation" && cat PROTEINS/all.annotation.tsv) \
    > temp &&
    mv temp PROTEINS/all.annotation.tsv

# check differences
cat PROTEINS/all.size.tsv |
    grep -F -f <(cut -f 1 PROTEINS/all.annotation.tsv) -v

tsv-join \
    PROTEINS/all.strain.tsv \
    --data-fields 1 \
    -f PROTEINS/all.size.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.strain_size.tsv

tsv-join \
    PROTEINS/all.strain_size.tsv \
    --data-fields 1 \
    -f PROTEINS/all.annotation.tsv \
    --key-fields 1 \
    --append-fields 2 \
    > PROTEINS/all.info.tsv

cat PROTEINS/all.info.tsv |
    wc -l
#7182873

```

## Phylogenetics with 40 single-copy genes

### Find corresponding proteins by `hmmsearch`

* Download HMM models as described in [`hmm/README.md`](../hmm/README.md)

* The `E_VALUE` was manually adjusted to 1e-20 to reach a balance between sensitivity and
  speciality.

```shell
E_VALUE=1e-20

cd ~/data/Pseudomonas

## example
#gzip -dcf ASSEMBLY/Ac_axa_ATCC_25176/*_protein.faa.gz |
#    hmmsearch -E 1e-20 --domE 1e-20 --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/BA00001.hmm - |
#    grep '>>' |
#    perl -nl -e '/>>\s+(\S+)/ and print $1'

# Find all genes
for marker in BA000{01..40}; do
    echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    for GENUS in $(cat genus.lst); do
        echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/scg40/bacteria_and_archaea_dir/${marker}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ and printf qq{%s\t%s\n}, $1, $ENV{STRAIN};
                '
        done \
            > PROTEINS/${marker}/${GENUS}.replace.tsv
    done

    echo
done


```

### Create a valid marker gene list

* `hmmsearch` may identify more than one copy for some marker genes.

    * BA00004: translation initiation factor EF-2
    * BA00005: translation initiation factor IF-2
    * BA00008: signal recognition particle protein

* Acinetobacter
    * BA00028

Compare proteins and strains.

```shell script
cd ~/data/Pseudomonas

for marker in BA000{01..03} BA000{06..07} BA000{09..27} BA000{29..40}; do
    echo ${marker}
done > marker.lst

for marker in $(cat marker.lst); do
    echo "==> marker [${marker}]"

    for GENUS in $(cat genus.lst); do
        cat PROTEINS/${marker}/${GENUS}.replace.tsv |
            cut -f 2 |
            diff - taxon/${GENUS}
    done

    echo
done

```

### Align and concat marker genes to create species tree

* Strains within a species share a large proportion of identical protein sequences.

* Use `trimal -automated1` to remove gaps.

```shell script
cd ~/data/Pseudomonas

# extract sequences
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for GENUS in $(cat genus.lst); do
            >&2 echo "==> GENUS [${GENUS}]"

            cat PROTEINS/{}/${GENUS}.replace.tsv |
                cut -f 1 |
                tsv-uniq |
                samtools faidx PROTEINS/all.uniq.fa -r -
        done \
            > PROTEINS/{}/{}.pro.fa

        for GENUS in $(cat genus.lst); do
            cat PROTEINS/{}/${GENUS}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        >&2 echo
    '

# aligning each markers with muscle
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> {}"

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '

for marker in $(cat marker.lst); do
    >&2 echo "==> marker [${marker}]"

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done

# concat marker genes
for marker in $(cat marker.lst); do
    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/scg40.aln.fas

fasops concat PROTEINS/scg40.aln.fas strains.lst -o PROTEINS/scg40.aln.fa

# trim with TrimAl
trimal -in PROTEINS/scg40.aln.fa -out PROTEINS/scg40.trim.fa -automated1

# FastTree produces NJ trees to simulate ML ones
FastTree PROTEINS/scg40.trim.fa > PROTEINS/scg40.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Pseudomonas/tree

nw_reroot ../PROTEINS/scg40.trim.newick B_sub_subtilis_168 St_aur_aureus_NCTC_8325 \
    > scg40.reroot.newick

rm scg40.condensed.map

# genus
bash ~/Scripts/withncbi/taxon/condense_tree.sh scg40.reroot.newick ../strains.taxon.tsv 1 4

mv scg40.condense.newick scg40.genus.newick
cat scg40.condense.map >> scg40.condensed.map

# species_group
bash ~/Scripts/withncbi/taxon/condense_tree.sh scg40.genus.newick ../strains.taxon.tsv 1 5

mv scg40.genus.condense.newick scg40.species_group.newick
cat scg40.genus.condense.map >> scg40.condensed.map

# species
bash ~/Scripts/withncbi/taxon/condense_tree.sh scg40.species_group.newick ../strains.taxon.tsv 1 6

mv scg40.species_group.condense.newick scg40.species.newick
cat scg40.species_group.condense.map >> scg40.condensed.map

rm *.condense.map

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 scg40.species.newick |
    rsvg-convert -o Pseudomonas.scg40.png

```

## Phylogenetics with bac120

### Find corresponding proteins by `hmmsearch`

```shell
E_VALUE=1e-20

cd ~/data/Pseudomonas

# Find all genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    echo "==> marker [${marker}]"

    mkdir -p PROTEINS/${marker}

    for GENUS in $(cat genus.lst); do
        echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/bac120/HMM/${marker}.HMM - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ and printf qq{%s\t%s\n}, $1, $ENV{STRAIN};
                '
        done \
            > PROTEINS/${marker}/${GENUS}.replace.tsv
    done

    echo
done


```

### Align and concat marker genes to create species tree

```shell script
cd ~/data/Pseudomonas

# extract sequences
cat ~/data/HMM/bac120/bac120.tsv |
    sed '1d' |
    cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for GENUS in $(cat genus.lst); do
            >&2 echo "==> GENUS [${GENUS}]"

            cat PROTEINS/{}/${GENUS}.replace.tsv |
                cut -f 1 |
                tsv-uniq |
                samtools faidx PROTEINS/all.uniq.fa -r -
        done \
            > PROTEINS/{}/{}.pro.fa

        for GENUS in $(cat genus.lst); do
            cat PROTEINS/{}/${GENUS}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        >&2 echo
    '

# aligning each markers with muscle
cat ~/data/HMM/bac120/bac120.tsv |
    sed '1d' |
    cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> {}"

        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    '

for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    >&2 echo "==> marker [${marker}]"

    # 1 name to many names
    cat PROTEINS/${marker}/${marker}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.aln.fa <(echo {}) stdout
        " \
        > PROTEINS/${marker}/${marker}.replace.fa
done

# concat marker genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/bac120.aln.fas

fasops concat PROTEINS/bac120.aln.fas strains.lst -o PROTEINS/bac120.aln.fa

# trim with TrimAl
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Pseudomonas/tree

nw_reroot ../PROTEINS/bac120.trim.newick B_sub_subtilis_168 St_aur_aureus_NCTC_8325 \
    > bac120.reroot.newick

rm bac120.condensed.map

# genus
bash ~/Scripts/withncbi/taxon/condense_tree.sh bac120.reroot.newick ../strains.taxon.tsv 1 4

mv bac120.condense.newick bac120.genus.newick
cat bac120.condense.map >> bac120.condensed.map

# species_group
bash ~/Scripts/withncbi/taxon/condense_tree.sh bac120.genus.newick ../strains.taxon.tsv 1 5

mv bac120.genus.condense.newick bac120.species_group.newick
cat bac120.genus.condense.map >> bac120.condensed.map

# species
bash ~/Scripts/withncbi/taxon/condense_tree.sh bac120.species_group.newick ../strains.taxon.tsv 1 6

mv bac120.species_group.condense.newick bac120.species.newick
cat bac120.species_group.condense.map >> bac120.condensed.map

rm *.condense.map

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 bac120.species.newick |
    rsvg-convert -o Pseudomonas.bac120.png

```
