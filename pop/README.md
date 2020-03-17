# Build alignments on an whole Eukaryotes genus

Or order, family, species.

Genus *Trichoderma* as example.

[TOC levels=1-3]: # ""

- [Build alignments on an whole Eukaryotes genus](#build-alignments-on-an-whole-eukaryotes-genus)
  - [Section 1: select strains and download sequences.](#section-1-select-strains-and-download-sequences)
    - [`pop/trichoderma.*.tsv`](#poptrichodermatsv)
    - [`wgs_prep.pl`](#wgs_preppl)
    - [`assembly_prep.pl`](#assembly_preppl)
  - [Section 2: create a raw phylogenetic tree by MinHash](#section-2-create-a-raw-phylogenetic-tree-by-minhash)
  - [Section 3: prepare sequences for `egaz`](#section-3-prepare-sequences-for-egaz)
    - [Manually](#manually)
    - [`egaz template --prep`](#egaz-template---prep)
  - [Section 4: generate alignments](#section-4-generate-alignments)
  - [Section 5: cleaning](#section-5-cleaning)
  - [FAQ](#faq)


## Section 1: select strains and download sequences.


Create `pop/trichoderma.wgs.tsv` manually. Names should only contain alphanumeric characters and
underscores. Be careful with tabs and spaces, because .tsv stands for Tab-separated values, white
spaces matters.

Check NCBI pages

* https://www.ncbi.nlm.nih.gov/Traces/wgs/?view=wgs&search=Trichoderma
* http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5543
* http://www.ncbi.nlm.nih.gov/genome/?term=txid5543[Organism:exp]
* http://www.ncbi.nlm.nih.gov/assembly?term=txid5543[Organism:exp]

And query a local `ar_genbank` DB. This is just a convenient but not accurate approach, especially
for sub-species parts.

```mysql
SELECT
    CONCAT(LEFT(genus, 1),
            LEFT(TRIM(REPLACE(species, genus, '')),
                3),
            REPLACE((REPLACE(organism_name, species, '')),
                ' ',
                '_')),
    SUBSTRING(wgs_master, 1, 4),
    organism_name,
    assembly_level
FROM
    ar_genbank.ar
WHERE
    genus = 'Trichoderma'
        AND wgs_master LIKE '%000%'
ORDER BY assembly_level, organism_name

```

For genus contains many species, you should be careful that "Gspe" (*G*enus *spe*cies) style
abbreviation may mix up two or more species.

```mysql
SELECT DISTINCT
    species
FROM
    ar_genbank.ar
WHERE
    genus = 'Trichoderma'

```

### `pop/trichoderma.*.tsv`

When the two approaches get very different number of strains, you run the following steps. Check
intermediate results on necessary.

Working directory should be `~/data/alignment/trichoderma` in this section.

```bash
export RANK_LEVEL=genus
export RANK_ID=5543
export RANK_NAME=trichoderma

mkdir -p ~/data/alignment/${RANK_NAME}            # Working directory
cd ~/data/alignment/${RANK_NAME}

```

#### `.wgs.tsv`

You can copy & paste the following block of codes.

```bash
# stage1
# Results from sql query.
mysql -ualignDB -palignDB gr_euk -e "
    SELECT 
        SUBSTRING(wgs,1,6) as prefix0,
        SUBSTRING(wgs,1,4) as prefix,
        organism_name,
        status 
    FROM gr 
    WHERE wgs != '' AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        CONCAT(SUBSTRING(wgs_master,1,5), RIGHT(wgs_master,1)) as prefix0,
        SUBSTRING(wgs_master,1,4) as prefix,
        organism_name,
        assembly_level 
    FROM ar 
    WHERE wgs_master != '' AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    >> raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT
        CONCAT(SUBSTRING(wgs_master,1,5), RIGHT(wgs_master,1)) as prefix0,
        SUBSTRING(wgs_master,1,4) as prefix,
        organism_name,
        assembly_level 
    FROM ar 
    WHERE wgs_master != '' AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    >> raw.tsv

# stage2
# NCBI changed its WGS page, make curl defunct
# Click on the 'Taxonomic Groups' on the left panel
# Click the 'Download' button in the middle of NCBI WGS page.
rm -f ~/Downloads/wgs_selector*.csv
chrome "https://www.ncbi.nlm.nih.gov/Traces/wgs/?view=wgs&search=${RANK_NAME}"

# Quit chrome Cmd-Q

# There're no chromosome level assemblies in WGS
cat ~/Downloads/wgs_selector.csv |
    perl -nl -a -F"," -e '
        my $p = substr($F[0],0,4);
        my $status = $F[13] > 0 ? q{Scaffold} : q{Contig};
        print qq{$F[0]\t$p\t$F[4]\t$status}
    ' \
    >> raw.tsv

cat raw.tsv |
    perl -nl -a -F"\t" -e '
        BEGIN{my %seen}; 
        /^prefix/i and next;
        scalar @F == 4 or next;
        $seen{$F[1]}++;
        $seen{$F[1]} > 1 and next;
        print join(qq{\t}, $F[0], $F[0], $F[2], $F[3]);
    ' \
    > raw2.tsv

# stage3
# Run `wgs_prep.pl` to get a crude `raw2.csv`
perl ~/Scripts/withncbi/taxon/wgs_prep.pl -f raw2.tsv --csvonly

echo -e '#name\tprefix\torganism\tcontigs' > raw3.tsv
cat raw2.csv |
    perl -nl -a -F"," -e '
        /^prefix/i and next;
        s/"//g for @F;
        @O = split(/ /, $F[3]);
        $name = substr($O[0],0,1) . substr($O[1],0,3);
        $name .= q{_} . $F[4] if $F[4];
        $name =~ s/\s+$//g;
        $name =~ s/\W+/_/g;
        print qq{$name\t$F[0]\t$F[3]\t$F[9]}
    ' |
    sort -t$'\t' -k4 -n |
    uniq \
    >> raw3.tsv

mv raw3.tsv ${RANK_NAME}.wgs.tsv

# find potential duplicated strains or assemblies
cat ${RANK_NAME}.wgs.tsv |
    perl -nl -a -F"\t" -e 'print $F[0]' |
    sort |
    uniq -c |
    sort -nr

# Edit .tsv, remove duplicated strains, check strain names and comment out poor assemblies.
# vim ${GENUS}.wgs.tsv

# Cleaning
rm raw*.*sv

```

#### `.assembly.tsv`

```bash

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
#        AND wgs_master = ''
#        AND assembly_level = 'Chromosome'
        AND organism_name != species
        AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT 
        organism_name, species, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
#        AND wgs_master = ''
#        AND assembly_level = 'Chromosome'
        AND organism_name != species
        AND ${RANK_LEVEL}_id = ${RANK_ID}
    " \
    >> raw.tsv

echo -e '#name\tftp_path\torganism\tassembly_level' > ${RANK_NAME}.assembly.tsv

cat raw.tsv |
    perl -nl -a -F"\t" -e '
        BEGIN{my %seen}; 
        /^organism_name/i and next;
        $n = $F[0];
        $rx = quotemeta $F[1];
        $n =~ s/$rx\s*//;
        $n =~ s/\s+$//;
        $n =~ s/\W+/_/g;
        @O = split(/ /, $F[1]);
        $name = substr($O[0],0,1) . substr($O[1],0,3);
        $name .= q{_} . $n if $n;
        $seen{$name}++;
        $seen{$name} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $name, $F[2], $F[1], $F[3];
        ' \
    >> ${RANK_NAME}.assembly.tsv


# comment out unneeded conditions

# find potential duplicated strains or assemblies
cat ${RANK_NAME}.assembly.tsv |
    cut -f 1 |
    sort |
    uniq -c |
    sort -nr

# Edit .tsv, remove unnecessary strains, check strain names and comment out poor assemblies.
# vim ${GENUS}.assembly.tsv

# Cleaning
rm raw*.*sv

```

```bash
unset RANK_LEVEL
unset RANK_ID
unset RANK_NAME

```

### `wgs_prep.pl`

Put the .tsv file to `~/Scripts/withncbi/pop/` and run `wgs_prep.pl` again. When everything is fine,
commit the .tsv file.

For detailed WGS info, click Prefix column lead to WGS project, where we could download gzipped
fasta and project description manually.

`wgs_prep.pl` will create a directory named `WGS` and some files containing meta information:

```bash
cd ~/data/alignment/trichoderma

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/withncbi/pop/trichoderma.wgs.tsv \
    --fix \
    -o WGS

```

1. `trichoderma.wgs.csv`

   Various information for WGS projects extracted from NCBI WGS record pages.

2. `trichoderma.wgs.rsync.sh` and `trichoderma.wgs.aria2.sh`

   Download WGS files.

3. `trichoderma.wgs.data.yml` and `trichoderma.wgs.master.yml`

   ```yaml
   ---
   data:
     - taxon: 452589
       name: 'Tart_IMI_2206040'
       sciname: 'Trichoderma atroviride IMI 206040'
       prefix: 'ABDG02'
       coverage: '8.26x Sanger'
   ```

Download WGS files.

```bash
cd ~/data/alignment/trichoderma

bash WGS/trichoderma.wgs.rsync.sh
bash WGS/trichoderma.wgs.aria2.sh

# check downloaded .gz files
find WGS -name "*.gz" | parallel -j 4 gzip -t

```

### `assembly_prep.pl`

Information of assemblies are collected from *_assembly_report.txt *after* downloading.

**Caution**: line endings of *_assembly_report.txt files are `CRLF`.

```bash
cd ~/data/alignment/trichoderma

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/trichoderma.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/trichoderma.assembly.rsync.sh

bash ASSEMBLY/trichoderma.assembly.collect.sh

```

## Section 2: create a raw phylogenetic tree by MinHash

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

