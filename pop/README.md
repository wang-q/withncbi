# Build alignments across a eukaryotic taxonomy rank

Genus *Trichoderma* as an example.

[TOC levels=1-3]: # ""

- [Build alignments on an whole Eukaryotes genus](#build-alignments-on-an-whole-eukaryotes-genus)
  - [Section 1: select strains and download sequences.](#section-1-select-strains-and-download-sequences)
    - [`pop/trichoderma.*.tsv`](#poptrichodermatsv)
    - [`wgs_prep.pl`](#wgs_preppl)
    - [`assembly_prep.pl`](#assembly_preppl)
  - [Section 2: create a raw phylogenetic tree by MinHash](#section-2-create-a-raw-phylogenetic-tree-by-minhash)
    - [`mash`](#mash)
  - [Section 3: prepare sequences for `egaz`](#section-3-prepare-sequences-for-egaz)
    - [Manually](#manually)
    - [`egaz template --prep`](#egaz-template---prep)
  - [Section 4: generate alignments](#section-4-generate-alignments)
  - [Section 5: cleaning](#section-5-cleaning)
  - [FAQ](#faq)


## Preparations

* Install `nwr` and create a local taxonomy database.

```shell
brew install wang-q/tap/nwr

nwr download
nwr txdb

```

* Create the [assembly database](https://github.com/wang-q/nwr/blob/master/doc/assembly.md).

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

Check also the outgroups of the family Hypocreaceae.

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


## Trichoderma: assembly

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

### `assembly_prep.pl`

Information of assemblies are collected from *_assembly_report.txt *after* downloading.

**Caution**: line endings of *_assembly_report.txt files are `CRLF`.

```bash
cd ~/data/alignment/trichoderma

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Trichoderma.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/trichoderma.assembly.rsync.sh

bash ASSEMBLY/trichoderma.assembly.collect.sh

```

## Section 2: create a raw phylogenetic tree by MinHash

### `mash`

```bash
mkdir -p ~/data/alignment/trichoderma/mash
cd ~/data/alignment/trichoderma/mash

for dir in $(find ../ASSEMBLY ../WGS -maxdepth 1 -mindepth 1 -type d ); do
    2>&1 echo "==> ${dir}"

    name=$(basename ${dir})

    if [[ -e ${name}.msh ]]; then
        continue
    fi

    find ${dir} -name "*.fsa_nt.gz" -or -name "*_genomic.fna.gz" |
        xargs cat |
        mash sketch -k 21 -s 100000 -p 4 - -I "${name}" -o ${name}
done

mash triangle -E -l <( find . -maxdepth 1 -type f -name "*.msh" | sort ) > dist.tsv

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

        group <- cutree(clusters, k=3) # h=0.5
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        cat(format_tsv(groups))
    '

nw_display -w 600 -b 'visibility:hidden' -s tree.nwk |
    rsvg-convert -o ~/Scripts/withncbi/image/Trichoderma.png

```

![Trichoderma.png](../image/Trichoderma.png)

## Section 3: prepare sequences for `egaz`

### Manually

* `perseq` mean split fasta by names, target or good assembles should set it
* `--species Fungi` specify the species or clade of this group for RepeatMasker

```bash
cd ~/data/alignment/trichoderma

mkdir -p GENOMES

# Sanger
for perseq in Tree_QM6a Tvir_Gv29_8 Tatr_IMI_206040; do
    echo ASSEMBLY/${perseq};
done |
    parallel --no-run-if-empty --linebuffer -k -j 2 '
        echo >&2 "==> {/}"

        if [ -d GENOMES/{/} ]; then
            echo >&2 "    GENOMES/{/} presents"
            exit;
        fi

        FILE_FA=$(ls {} | grep "_genomic.fna.gz" | grep -v "_from_")
        echo >&2 "==> Processing {}/${FILE_FA}"

        egaz prepseq \
            {}/${FILE_FA} \
            -o GENOMES/{/} \
            --min 50000 --gi -v --repeatmasker "--species Fungi --parallel 8"

        FILE_GFF=$(ls {} | grep "_genomic.gff.gz")
        echo >&2 "==> Processing {}/${FILE_GFF}"

        gzip -d -c {}/${FILE_GFF} > GENOMES/{/}/chr.gff
    '

# Other assemblies
find ASSEMBLY -maxdepth 1 -type d -path "*/????*" |
    parallel --no-run-if-empty --linebuffer -k -j 2 '
        echo >&2 "==> {/}"

        if [ -d GENOMES/{/} ]; then
            echo >&2 "    GENOMES/{/} presents"
            exit;
        fi

        FILE_FA=$(ls {} | grep "_genomic.fna.gz" | grep -v "_from_")

        egaz prepseq \
            {}/${FILE_FA} \
            -o GENOMES/{/} \
            --about 5000000 \
            --min 5000 --gi -v --repeatmasker "--species Fungi --parallel 8"

        FILE_GFF=$(ls {} | grep "_genomic.gff.gz")
        echo >&2 "==> Processing {}/${FILE_GFF}"

        gzip -d -c {}/${FILE_GFF} > GENOMES/{/}/chr.gff
    '

# WGS
find WGS -maxdepth 1 -type d -path "*/????*" |
    parallel --no-run-if-empty --linebuffer -k -j 2 '
        echo >&2 "==> {/}"

        if [ -d GENOMES/{/} ]; then
            echo >&2 "    GENOMES/{/} presents"
            exit;
        fi

        FILE_FA=$(ls {} | grep ".fsa_nt.gz")

        egaz prepseq \
            {}/${FILE_FA} \
            -o GENOMES/{/} \
            --about 5000000 \
            --min 5000 --gi -v --repeatmasker "--species Fungi --parallel 8"
    '


```

### `egaz template --prep`

Or use `egaz template --prep`. In this approach, GFF files should be manually placed in the GENOMES/
directory.

```bash
cd ~/data/alignment/trichoderma

egaz template \
    ASSEMBLY WGS \
    --prep -o GENOMES \
    --perseq Tree_QM6a --perseq Tvir_Gv29_8 --perseq Tatr_IMI_206040 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Fungi --parallel 8"

bash GENOMES/0_prep.sh

```

## Section 4: generate alignments

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

## Section 5: cleaning

This is the ultimate final step. No more. Actually you may choose not to do this. It's depended on
your disk capacity.

```bash
cd ~/data/alignment/Fungi/trichoderma

# clean raw fasta
find . -maxdepth 1 -type d -name "*_raw" | xargs rm -fr

# clean maf-fasta
find . -maxdepth 1 -type d -name "*_fasta" | xargs rm -fr

# clean raxml phy
find . -maxdepth 2 -type f -name "*.phy" -or -name "*.phy.reduced" | xargs rm

# compress files
find . -type f -name "*.maf" | parallel gzip
find . -type f -name "*.fas" | parallel gzip

```

## FAQ

* Why .tsv? All of your other programs use .csv.

  There are strains of which sub-species parts contain commas, can you believe it?

* I've 500 genomes of *E. coli*, the manually editing step 1 kills me.

  The whole `pop/` and much of `util/` scripts are for Eukaryotes. For small genomes of bacteria,
  archaea and organelles, check `taxon/bacteria_gr.md`.

* Your command lines executed and the results are wired.

  Be sure you aren't in Windows. Be sure you are familiar to bash command lines.

  Or send what you want to me and let me do the job.

* I have a very good assembly on chromosome level, but I can't find it in WGS.

  Best genomes on the world went to NCBI RefSeq. Use tools in `util/` to download them. Examples can
  be found in `pop/OPs-download.md`.

