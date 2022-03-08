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
    * [Phylogenetics with 40 single-copy genes](#phylogenetics-with-40-single-copy-genes)
        + [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
        + [Create a valid marker gene list](#create-a-valid-marker-gene-list)
        + [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
        + [Tweak the concat tree](#tweak-the-concat-tree)
    * [Phylogenetics with bac120](#phylogenetics-with-bac120)
        + [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch-1)
        + [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree-1)
        + [Tweak the concat tree](#tweak-the-concat-tree-1)
    * [Protein domains and families](#protein-domains-and-families)
        + [Proteins in Pseudomonas strains](#proteins-in-pseudomonas-strains)
        + [Scrap PFAM domains](#scrap-pfam-domains)
        + [Scan every domain](#scan-every-domain)
        + [InterProScan](#interproscan)

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
brew install brewsci/bio/muscle
brew install brewsci/bio/fasttree
brew install brewsci/bio/easel
brew install brewsci/bio/newick-utils
brew install brewsci/bio/trimal

brew install datamash
brew install miller
brew install wang-q/tap/tsv-utils

brew install librsvg
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

# Pseudomonas aeruginosa PAO1 is in the reference list
cat reference.tsv |
    tsv-select -H -f organism_name,species,genus,ftp_path,assembly_level \
    > raw.tsv

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
        $seen{$F[3]}++; # ftp_path
        $seen{$F[3]} > 1 and next;
        $seen{$F[5]}++; # abbr_name
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

* Some strains were anomalously labeled and identified by the `mash` tree.
    * Pseudom_flu_GCF_900636635_1
    * Pseudom_chl_GCF_001023535_1
    * Pseudom_syr_GCF_004006335_1
    * Pseudom_puti_GCF_003228315_1 and Pseudom_puti_GCF_020172705_1

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
#  1520 ASSEMBLY/n50.pass.csv
#  1520 ASSEMBLY/n50.tsv

tsv-join \
    ASSEMBLY/Pseudomonas.assembly.collect.csv \
    --delimiter "," -H --key-fields 1 \
    --filter-file ASSEMBLY/n50.pass.csv |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_900636635 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_001023535 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_004006335 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_003228315 |
    tsv-filter --delimiter "," -H --str-not-in-fld 1:GCF_020172705 \
    > ASSEMBLY/Pseudomonas.assembly.pass.csv

wc -l ASSEMBLY/Pseudomonas.assembly*csv
#  1520 ASSEMBLY/Pseudomonas.assembly.collect.csv
#  1515 ASSEMBLY/Pseudomonas.assembly.pass.csv

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
| 469     | Acinetobacter     | 43       | 483      |
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
| 286     | Pseudomonas       | 172      | 824      |
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

nw_reroot ../mash/tree.nwk B_sub_subtilis_168 St_aur_aureus_NCTC_8325 |
    nw_order -c n - \
    > mash.reroot.newick

# rank::col
ARRAY=(
#    'order::7'
    'family::6'
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
    rsvg-convert -o Pseudomonas.mash.png

```

## Collect proteins

### `all.pro.fa`

```shell script
cd ~/data/Pseudomonas

mkdir -p PROTEINS

find ASSEMBLY -maxdepth 1 -mindepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l
# 1519

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l
# 1519

cat strains.lst |
    wc -l
# 1514

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
#7151013

cat PROTEINS/all.pro.fa |
    grep "^>" |
    tsv-uniq |
    wc -l
#2468817

# annotations may be different
cat PROTEINS/all.uniq.fa |
    grep "^>" |
    wc -l
#2420143

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
#7151013

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
#7151013

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
#7151014

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

# Extract sequences
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for GENUS in $(cat genus.lst); do
            cat PROTEINS/{}/${GENUS}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with `muscle`
cat marker.lst |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

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

# Concat marker genes
for marker in $(cat marker.lst); do
    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/scg40.aln.fas

fasops concat PROTEINS/scg40.aln.fas strains.lst -o PROTEINS/scg40.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/scg40.aln.fa -out PROTEINS/scg40.trim.fa -automated1

# FastTree produces NJ trees to simulate ML ones
FastTree PROTEINS/scg40.trim.fa > PROTEINS/scg40.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Pseudomonas/tree

nw_reroot ../PROTEINS/scg40.trim.newick B_sub_subtilis_168 St_aur_aureus_NCTC_8325 |
    nw_order -c n - \
    > scg40.reroot.newick

# rank::col
ARRAY=(
#    'order::7'
#    'family::6'
    'genus::5'
    'species::4'
)

rm scg40.condensed.map
CUR_TREE=scg40.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick scg40.${GROUP_NAME}.newick
    cat condense.map >> scg40.condensed.map

    CUR_TREE=scg40.${GROUP_NAME}.newick
done

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

# Extract sequences
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

        for GENUS in $(cat genus.lst); do
            cat PROTEINS/{}/${GENUS}.replace.tsv
        done \
            > PROTEINS/{}/{}.replace.tsv

        faops some PROTEINS/all.uniq.fa <(
            cat PROTEINS/{}/{}.replace.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout \
            > PROTEINS/{}/{}.pro.fa
    '

# Align each markers with muscle
cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1 |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo "==> marker [{}]"

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

# Concat marker genes
for marker in $(cat ~/data/HMM/bac120/bac120.tsv | sed '1d' | cut -f 1); do
    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.replace.fa stdout

    # empty line for .fas
    echo
done \
    > PROTEINS/bac120.aln.fas

fasops concat PROTEINS/bac120.aln.fas strains.lst -o PROTEINS/bac120.aln.fa

# Trim poorly aligned regions with `TrimAl`
trimal -in PROTEINS/bac120.aln.fa -out PROTEINS/bac120.trim.fa -automated1

# To make it faster
FastTree -fastest -noml PROTEINS/bac120.trim.fa > PROTEINS/bac120.trim.newick

```

### Tweak the concat tree

```shell script
cd ~/data/Pseudomonas/tree

nw_reroot ../PROTEINS/bac120.trim.newick B_sub_subtilis_168 St_aur_aureus_NCTC_8325 |
    nw_order -c n - \
    > bac120.reroot.newick

rm bac120.condensed.map

# rank::col
ARRAY=(
#    'order::7'
#    'family::6'
    'genus::5'
    'species::4'
)

rm bac120.condensed.map
CUR_TREE=bac120.reroot.newick

for item in "${ARRAY[@]}" ; do
    GROUP_NAME="${item%%::*}"
    GROUP_COL="${item##*::}"

    bash ~/Scripts/withncbi/taxon/condense_tree.sh ${CUR_TREE} ../strains.taxon.tsv 1 ${GROUP_COL}

    mv condense.newick bac120.${GROUP_NAME}.newick
    cat condense.map >> bac120.condensed.map

    CUR_TREE=bac120.${GROUP_NAME}.newick
done

# png
nw_display -s -b 'visibility:hidden' -w 600 -v 30 bac120.species.newick |
    rsvg-convert -o Pseudomonas.bac120.png

```

## Protein domains and families

### Proteins in Pseudomonas strains

* [`GO:0005975` carbohydrate metabolic process](https://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0005975)

* Search `https://www.pseudomonas.com/goterms` for `GO:0005975`
    * 107 - Pseudomonas aeruginosa PAO1
    * 7360 - Pseudomonas putida KT2440
    * 479 - Pseudomonas chlororaphis subsp. aureofaciens 30-84
    * 116 - Pseudomonas fluorescens SBW25
    * 113 - Pseudomonas protegens Pf-5
    * 123 - Pseudomonas stutzeri A1501
    * 112 - Pseudomonas syringae pv. syringae B728a
    * 114 - Pseudomonas savastanoi pv. phaseolicola 1448A
    * 117 - Pseudomonas entomophila L48
    * 109 - Pseudomonas aeruginosa UCBPP-PA14

```shell
cd ~/data/Pseudomonas/
mkdir -p ~/data/Pseudomonas/DOMAINS/ref_strain

for ID in 107 7360 479 116 113 123 112 114 117 109 ; do
    URL=$(printf 'https://www.pseudomonas.com/goterms/list?accession=GO:0005975&strain_id=%d&format=TAB' $ID)
    curl -L ${URL}
done  |
    sed '1d' |
    tsv-select -f 1 \
    > DOMAINS/ref_strain/locus.lst

gzip -dcf \
    ASSEMBLY/Pseudom_aer_PAO1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_puti_KT2440_GCF_000007565_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_chl_aureofaciens_30_84_GCF_000281915_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_flu_SBW25_GCF_000009225_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_pro_Pf_5_GCF_000012265_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_stu_A1501_GCF_000013785_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_syr_pv_syringae_B728a_GCF_000012245_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_ento_L48_GCF_000026105_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_aer_UCBPP_PA14_GCF_000014625_1/*_genomic.gff.gz |
    grep -v "^#" |
    grep -F -w -f DOMAINS/ref_strain/locus.lst |
    tsv-filter --str-eq 3:gene |
    perl -nl -e 'print $1 if /\bID=(.+?);/i' \
    > DOMAINS/ref_strain/gene.lst

gzip -dcf \
    ASSEMBLY/Pseudom_aer_PAO1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_puti_KT2440_GCF_000007565_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_chl_aureofaciens_30_84_GCF_000281915_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_flu_SBW25_GCF_000009225_2/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_pro_Pf_5_GCF_000012265_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_stu_A1501_GCF_000013785_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_syr_pv_syringae_B728a_GCF_000012245_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_ento_L48_GCF_000026105_1/*_genomic.gff.gz \
    ASSEMBLY/Pseudom_aer_UCBPP_PA14_GCF_000014625_1/*_genomic.gff.gz |
    grep -v "^#" |
    grep -F -w -f DOMAINS/ref_strain/gene.lst |
    tsv-filter --str-eq 3:CDS |
    perl -nl -e 'print $1 if /\bName=(.+?);/i' \
    > DOMAINS/ref_strain/pro.lst

wc -l DOMAINS/ref_strain/*.lst
#  497 DOMAINS/ref_strain/gene.lst
#  513 DOMAINS/ref_strain/locus.lst
#  497 DOMAINS/ref_strain/pro.lst

faops some PROTEINS/all.uniq.fa DOMAINS/ref_strain/pro.lst stdout \
    > DOMAINS/ref_strain/pro.fa

faops size DOMAINS/ref_strain/pro.fa | wc -l
# 495

# prepare an HMM database for faster hmmscan searches
#hmmpress ~/data/HMM/PFAM/Pfam-A.hmm

E_VALUE=1e-5

hmmscan --cpu 4 -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw \
    ~/data/HMM/PFAM/Pfam-A.hmm DOMAINS/ref_strain/pro.fa |
    grep '>>' |
    perl -nl -e '/>>\s+(\S+)/ and print $1' |
    tee DOMAINS/ref_strain/domain.lst

# query ID against .hmm.dat to find AC
cat ~/data/HMM/PFAM/Pfam-A.hmm.dat |
    grep -F -w -f DOMAINS/ref_strain/domain.lst -A 2 |
    grep -E ' (ID|AC) ' |
    perl -nl -e 'print substr($_, 10) ' |
    paste -d $'\t' - - |
    perl -nlp -e 's/\.\d+$//g' |
    tsv-select -f 2,1 \
    > DOMAINS/ref_strain/pfam_domain.tsv

wc -l < DOMAINS/ref_strain/pfam_domain.tsv
# 107


```

### Scrap PFAM domains

* Perform keyword search in `pfam` and save the result page as `html only`.
    * [GO:0005975](https://pfam.xfam.org/search/keyword?query=GO%3A0005975)
    * [Glyco_hyd](http://pfam.xfam.org/search/keyword?query=Glyco_hyd)

```shell
cd ~/data/Pseudomonas/

cp DOMAINS/ref_strain/pfam_domain.tsv raw.tsv

cat GO_0005975.htm |
    pup 'table.resultTable tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - |
    tsv-select -f 2-4 \
    >> raw.tsv

cat Glyco_hyd.htm |
    pup 'table.resultTable tr td text{}' |
    grep '\S' |
    paste -d $'\t' - - - - - |
    tsv-select -f 2-4 \
    >> raw.tsv

#* [carbohydrate metabolic](http://pfam.xfam.org/search/keyword?query=carbohydrate+metabolic)
#cat carbohydrate_metabolic.htm |
#    pup 'table.resultTable tr td text{}' |
#    grep '\S' |
#    paste -d $'\t' - - - - - - - - |
#    tsv-select -f 2-4 \
#    >> raw.tsv

#* [Glycosyl_hydrolase](https://pfam.xfam.org/search/keyword?query=Glycosyl+hydrolase)
#cat Glycosyl_hydrolase.htm |
#    pup 'table.resultTable tr td text{}' |
#    grep '\S' |
#    paste -d $'\t' - - - - - - - - - |
#    tsv-filter --ne 5:10000000 | # Text fields of Pfam entries are not empty
#    tsv-select -f 2-4 \
#    >> raw.tsv

wc -l < raw.tsv
#290

cat raw.tsv |
    tsv-filter --str-not-in-fld 2:"DUF" |
    tsv-uniq -f 1 |
    tsv-sort -k2,2 \
    > pfam_domain.tsv

wc -l < pfam_domain.tsv
#216

mkdir -p ~/data/Pseudomonas/DOMAINS/HMM

cat pfam_domain.tsv |
    parallel --col-sep "\t" --no-run-if-empty --linebuffer -k -j 4 '
        >&2 echo {2}
        curl -L http://pfam.xfam.org/family/{1}/hmm > DOMAINS/HMM/{2}.hmm
    '

find DOMAINS/HMM -type f -name "*.hmm" |
    wc -l
#216

```

### Scan every domain

* The `E_VALUE` was adjusted to 1e-3 to capture all possible sequences.

```shell script
E_VALUE=1e-3

cd ~/data/Pseudomonas/

for domain in $(cat pfam_domain.tsv | cut -f 2 | sort); do
    >&2 echo "==> domain [${domain}]"

    if [ -e DOMAINS/${domain}.replace.tsv ]; then
        continue;
    fi

    for GENUS in $(cat genus.lst); do
        >&2 echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw DOMAINS/HMM/${domain}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done
    done \
        > DOMAINS/${domain}.replace.tsv

    >&2 echo
done

for domain in $(cat pfam_domain.tsv | cut -f 2 | sort); do
    wc -l DOMAINS/${domain}.replace.tsv
done |
    datamash reverse -W |
    tsv-filter --ge 2:2000 |
    tsv-sort -k2,2nr |
    (echo -e "Domain\tCount" && cat) |
    mlr --itsv --omd cat

```

| Domain                              | Count |
|-------------------------------------|-------|
| DOMAINS/Epimerase.replace.tsv       | 53528 |
| DOMAINS/AAA_33.replace.tsv          | 36056 |
| DOMAINS/Hydrolase.replace.tsv       | 30648 |
| DOMAINS/HAD.replace.tsv             | 24827 |
| DOMAINS/F420_oxidored.replace.tsv   | 23957 |
| DOMAINS/RmlD_sub_bind.replace.tsv   | 22072 |
| DOMAINS/GDP_Man_Dehyd.replace.tsv   | 20434 |
| DOMAINS/HAD_2.replace.tsv           | 19649 |
| DOMAINS/ApbA.replace.tsv            | 17803 |
| DOMAINS/CBS.replace.tsv             | 17485 |
| DOMAINS/Hydrolase_3.replace.tsv     | 16587 |
| DOMAINS/3Beta_HSD.replace.tsv       | 15744 |
| DOMAINS/Glycos_transf_2.replace.tsv | 14495 |
| DOMAINS/Glyco_tranf_2_3.replace.tsv | 13715 |
| DOMAINS/Hydrolase_like.replace.tsv  | 12776 |
| DOMAINS/Glyco_trans_4_4.replace.tsv | 11103 |
| DOMAINS/NAD_Gly3P_dh_N.replace.tsv  | 8591  |
| DOMAINS/PfkB.replace.tsv            | 8207  |
| DOMAINS/SIS.replace.tsv             | 8120  |
| DOMAINS/Glyco_trans_2_3.replace.tsv | 7641  |
| DOMAINS/AP_endonuc_2.replace.tsv    | 6976  |
| DOMAINS/Alpha-amylase.replace.tsv   | 6546  |
| DOMAINS/Phos_pyr_kin.replace.tsv    | 6081  |
| DOMAINS/FGGY_C.replace.tsv          | 5874  |
| DOMAINS/CTP_transf_like.replace.tsv | 5325  |
| DOMAINS/Polysacc_deac_1.replace.tsv | 5282  |
| DOMAINS/Glyco_transf_21.replace.tsv | 4734  |
| DOMAINS/SKI.replace.tsv             | 4232  |
| DOMAINS/PGM_PMM_I.replace.tsv       | 4120  |
| DOMAINS/PGM_PMM_II.replace.tsv      | 4010  |
| DOMAINS/PGM_PMM_III.replace.tsv     | 4008  |
| DOMAINS/PGM_PMM_IV.replace.tsv      | 4006  |
| DOMAINS/FGGY_N.replace.tsv          | 3754  |
| DOMAINS/QRPTase_C.replace.tsv       | 3668  |
| DOMAINS/DctQ.replace.tsv            | 3063  |
| DOMAINS/Aldose_epim.replace.tsv     | 2768  |
| DOMAINS/CBM_48.replace.tsv          | 2763  |
| DOMAINS/Glyco_hydro_3.replace.tsv   | 2743  |
| DOMAINS/LamB_YcsF.replace.tsv       | 2467  |
| DOMAINS/Glyco_transf_28.replace.tsv | 2147  |
| DOMAINS/Glyco_hydro_19.replace.tsv  | 2122  |
| DOMAINS/Chitin_synth_2.replace.tsv  | 2047  |

Check each domain, such as `https://pfam.xfam.org/family/Epimerase`. Some domains are not directly
related to carbohydrate metabolism.

* ATP, AMP, GDP, CTP, and NADP associated
    * AAA_33
    * GDP_Man_Dehyd
    * CBS*
    * CTP_transf_like
* Binding to other small molecules
    * HAD* 卤素
    * F420*
    * SKI 莽草酸
    * QRPTase_C 喹啉酸
* Oxidoreductases
    * Epimerase
    * ApbA
    * 3Beta_HSD
    * NAD_Gly3P_dh_N

* We can't interproscan all `2420143` proteins

```shell script
cd ~/data/Pseudomonas

# All proteins appeared
cat pfam_domain.tsv | cut -f 2 | sort |
    tsv-filter --not-regex 1:'^AAA' |
    tsv-filter --not-regex 1:'^GDP' |
    tsv-filter --not-regex 1:'^CBS' |
    tsv-filter --not-regex 1:'^CTP' |
    tsv-filter --not-regex 1:'^HAD' |
    tsv-filter --not-regex 1:'^F420' |
    tsv-filter --not-regex 1:'^SKI' |
    tsv-filter --not-regex 1:'^QRPTase' |
    tsv-filter --not-regex 1:'^Epimerase' |
    tsv-filter --not-regex 1:'^ApbA' |
    tsv-filter --not-regex 1:'^3Beta_HSD' |
    tsv-filter --not-regex 1:'^NAD_Gly3P' \
    > domain.lst

wc -l < domain.lst
#202

cat domain.lst |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        tsv-select -f 2 DOMAINS/{}.replace.tsv
    ' |
    sort -u \
    > DOMAINS/domains.tsv

wc -l < DOMAINS/domains.tsv
#168791

faops size PROTEINS/all.uniq.fa | wc -l
#2420143

for domain in $(cat domain.lst); do
    echo 1>&2 "==> domain [${domain}]"

    tsv-join \
        DOMAINS/domains.tsv \
        --data-fields 1 \
        -f <(
            cat DOMAINS/${domain}.replace.tsv |
                perl -nla -e 'print qq{$F[1]\tO}'
        ) \
        --key-fields 1 \
        --append-fields 2 \
        --write-all "" \
        > DOMAINS/tmp.tsv

    mv DOMAINS/tmp.tsv DOMAINS/domains.tsv
done

datamash check < DOMAINS/domains.tsv
#168791 lines, 203 fields

# Add header line
for domain in $(cat domain.lst); do
    echo "${domain}"
done |
    (echo -e "#name" && cat) |
    paste -s -d $'\t' - \
    > DOMAINS/header.tsv

cat DOMAINS/header.tsv DOMAINS/domains.tsv \
    > tmp.tsv && mv tmp.tsv DOMAINS/domains.tsv

tsv-join \
    PROTEINS/all.info.tsv \
    --data-fields 1 \
    -f DOMAINS/domains.tsv \
    --key-fields 1 \
    --append-fields 2-203 |
     keep-header -- sort -k1,1 \
    > tmp.tsv && mv tmp.tsv DOMAINS/domains.tsv

datamash check < DOMAINS/domains.tsv
#168792 lines, 206 fields

rm DOMAINS/header.tsv

```

### InterProScan

InterProScan 启动很慢, 因此一次提交整个菌株里的数十个蛋白.

```shell script
cd ~/data/Pseudomonas

mkdir -p IPS

# extract wanted sequences
for GENUS in $(cat genus.lst); do
    >&2 echo "==> GENUS [${GENUS}]"

    cat taxon/${GENUS} |
        parallel --no-run-if-empty --linebuffer -k -j 8 '
            mkdir -p IPS/{}

            cat PROTEINS/all.info.tsv |
                tsv-filter --str-eq 2:{} |
                cut -f 1 |
                grep -Fx -f <(cut -f 1 DOMAINS/domains.tsv | grep "^{}") \
                > IPS/{}/wanted.lst

            faops some PROTEINS/all.replace.fa IPS/{}/wanted.lst IPS/{}/{}.fa
        '
done

# scan proteins of each strain with InterProScan
# By default InterProScan uses 8 cpu cores
mkdir -p split
split -a 4 -l 30 -d strains.lst split/

for f in $(find split -maxdepth 1 -type f -name "[0-9]*" | sort); do
    >&2 echo "==> IPS [${f}]"
    bsub -q mpi -n 24 -J "IPS-${f}" "
        cat ${f} |
            parallel --no-run-if-empty --linebuffer -k -j 6 '
                if [[ -e IPS/{}/{}.tsv ]]; then
                    >&2 echo {};
                    exit;
                fi

                interproscan.sh --cpu 4 -dp -f tsv,json,svg -i IPS/{}/{}.fa --output-file-base IPS/{}/{}
            '
        "
done

rm -fr split output.*

# IPS family
for GENUS in $(cat genus.lst); do
    echo 1>&2 "==> GENUS [${GENUS}]"

    for STRAIN in $(cat taxon/${GENUS}); do
        cat IPS/${STRAIN}/${STRAIN}.json |
            jq .results |
            jq -r -c '
                .[] |
                .xref as $name |
                .matches[] |
                .signature.entry |
                select(.type == "FAMILY") |
                [$name[0].name, .accession, .description] |
                @tsv
            ' |
            tsv-uniq -f 1
    done
done |
    (echo -e "#name\tfamily\tdescription" && cat) \
    > IPS/family.tsv

tsv-join \
    <(cut -f 1-4 DOMAINS/domains.tsv) \
    --data-fields 1 \
    -f IPS/family.tsv \
    --key-fields 1 \
    --append-fields 2-3 \
    --write-all "" |
    tsv-join \
        --data-fields 1 \
        -f DOMAINS/domains.tsv \
        --key-fields 1 \
        --append-fields 5-206 |
     keep-header -- sort -k1,1 \
    > IPS/predicts.tsv

```

