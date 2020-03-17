# JGI data

Login [phytozome portal](https://phytozome.jgi.doe.gov/pz/portal.html).

Select needed species in *PhytozomeV12*
[download page](http://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=PhytozomeV12),
then click the `Download Selected Files` button. A file named `PhytozomeV12_download.zip` would be
created on the fly. Extract compressed file to `~/data/PhytozomeV12/SpeciesName`.

JGI recommends using Globus, but it's blocked by GFW.

Check identities between Ensembl and JGI assemblies. Then use gff3 files from JGI to replace ones
from Ensembl.

## Organism

| Organism                      | Version   | JGI abbrev.       | Taxon ID | Restricted |
|:------------------------------|:----------|:------------------|:---------|:----------:|
| Arabidopsis thaliana columbia | TAIR10    | Athaliana         | 3702     |            |
| Arabidopsis lyrata            | v2.1      | Alyrata           | 81972    |            |
| Arabidopsis halleri           | v1.1      | Ahalleri          | 81970    |     o      |
| Capsella grandiflora          | v1.1      | Cgrandiflora      | 264402   |            |
| Capsella rubella              | v1.0      | Crubella          | 81985    |            |
| Boechera stricta              | v1.2      | Bstricta          | 72658    |            |
| Brassica rapa FPsc            | v1.3      | BrapaFPsc         | 3711     |     o      |
| Brassica oleracea capitata    | v1.0      | Boleraceacapitata | 3716     |            |
| Eutrema salsugineum           | v1.0      | Esalsugineum      | 72664    |            |
| Carica papaya                 | ASGPBv0.4 | Cpapaya           | 3649     |            |
| Gossypium raimondii           | v2.1      | Graimondii        | 29730    |            |
| Theobroma cacao               | v1.1      | Tcacao            | 3641     |            |
| Oryza sativa                  | v7.0      | Osativa           | 39947    |            |
| Oryza sativa Kitaake          | v3.1 (ER) | OsativaKitaake    |          |            |
| Brachypodium stacei           | v1.1      | Bstacei           | 1071399  |            |
| Brachypodium distachyon       | v3.1      | Bdistachyon       | 15368    |            |

## Stats

```bash

for name in \
    Athaliana Alyrata Ahalleri \
    Cgrandiflora Crubella Bstricta \
    BrapaFPsc Boleraceacapitata \
    Esalsugineum Cpapaya \
    Graimondii Tcacao \
    Osativa  OsativaKitaake \
    Bstacei Bdistachyon \
    ; do
    1>&2 echo "==> ${name}"
    
    find ~/data/PhytozomeV12/${name}/assembly -type f -name "*.fa.gz" -not -name "*masked*" |
        xargs cat |
        faops n50 -C -S stdin |
        (echo -e "Name\t${name}" && cat) |
        datamash transpose
done |
    tsv-uniq |
    mlr --itsv --omd cat

```

| Name              | N50      | S         | C     |
|:------------------|:---------|:----------|:------|
| Athaliana         | 23459830 | 119667750 | 7     |
| Alyrata           | 24464547 | 206667935 | 695   |
| Ahalleri          | 29271    | 127615339 | 11241 |
| Cgrandiflora      | 112041   | 105346052 | 4997  |
| Crubella          | 15060676 | 134834574 | 853   |
| Bstricta          | 2187891  | 189344188 | 1990  |
| BrapaFPsc         | 28488603 | 315053614 | 5713  |
| Boleraceacapitata | 40895475 | 385006588 | 9     |
| Esalsugineum      | 13441892 | 243117811 | 639   |
| Cpapaya           | 1182040  | 342680090 | 5901  |
| Graimondii        | 62175169 | 761406121 | 1033  |
| Tcacao            | 34397752 | 346164918 | 713   |
| Osativa           | 29958434 | 374471240 | 14    |
| OsativaKitaake    | 30273398 | 381570803 | 33    |
| Bstacei           | 23060899 | 234142426 | 112   |
| Bdistachyon       | 59130575 | 271163419 | 10    |

## MinHash

```bash
mkdir -p ~/data/alignment/JGI/mash
cd ~/data/alignment/JGI/mash

for name in \
    Athaliana Alyrata Ahalleri \
    Cgrandiflora Crubella Bstricta \
    BrapaFPsc Boleraceacapitata \
    Esalsugineum Cpapaya \
    Graimondii Tcacao \
    Osativa  OsativaKitaake \
    Bstacei Bdistachyon \
    ; do
    1>&2 echo "==> ${name}"
    
    find ~/data/PhytozomeV12/${name}/assembly -type f -name "*.fa.gz" -not -name "*masked*" |
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
        
        group <- cutree(clusters, k=2) # h=0.2
        groups <- as.data.frame(group)
        groups$ids <- rownames(groups)
        rownames(groups) <- NULL
        groups <- groups[order(groups$group), ]
        cat(format_tsv(groups))
    '

```

## prepseq

Skip Ahalleri, Cgrandiflora.

```bash
mkdir -p ~/data/alignment/JGI
cd ~/data/alignment/JGI

# sequences
for name in \
    Athaliana Alyrata \
    Crubella Bstricta \
    BrapaFPsc Boleraceacapitata \
    Esalsugineum Cpapaya \
    Graimondii Tcacao \
    Osativa  OsativaKitaake \
    Bstacei Bdistachyon \
    ; do
    
    if [ ! -d ~/data/alignment/JGI/${name} ]; then
        1>&2 echo "==> ${name}"
    
        mkdir -p ~/data/alignment/JGI/${name}
        pushd ~/data/alignment/JGI/${name} > /dev/null
        
        find ~/data/PhytozomeV12/${name}/assembly -type f -name "*.fa.gz" -not -name "*masked*" |
            xargs gzip -dc > toplevel.fa
        
        faops count toplevel.fa |
            perl -nla -e '
                next if $F[0] eq 'total';
                print $F[0] if $F[1] > 1000000;
                print $F[0] if $F[1] > 100000 and $F[6]/$F[1] < 0.05;
            ' |
            uniq > listFile
        faops some toplevel.fa listFile toplevel.filtered.fa
        faops split-name toplevel.filtered.fa .
        rm toplevel.fa toplevel.filtered.fa listFile

        popd  > /dev/null
    else
        1>&2 echo "==> ~/data/alignment/JGI/${name} exists"
    fi

done

# RepeatMasker
rsync -avP ~/data/alignment/JGI/ wangq@202.119.37.251:data/alignment/JGI

for name in \
    Athaliana Alyrata \
    Crubella Bstricta \
    BrapaFPsc Boleraceacapitata \
    Esalsugineum Cpapaya \
    Graimondii Tcacao \
    Osativa  OsativaKitaake \
    Bstacei Bdistachyon \
    ; do
    
    if [ ! -f ~/data/alignment/JGI/${name}/chr.2bit ]; then
        1>&2 echo "==> ${name}"
    
#        egaz prepseq \
#            --repeatmasker '--species Viridiplantae --gff --parallel 8' -v \
#            ~/data/alignment/JGI/${name}
        bsub -q mpi -n 24 -J "prep-${name}" \
            "egaz prepseq --repeatmasker '--species Viridiplantae --gff --parallel 24' -v ~/data/alignment/JGI/${name}"
    fi

done

# bjobs -w
# bkill -J "prep-*"

rsync -avP wangq@202.119.37.251:data/alignment/JGI/ ~/data/alignment/JGI

# gff
for name in \
    Athaliana Alyrata \
    Crubella Bstricta \
    BrapaFPsc Boleraceacapitata \
    Esalsugineum Cpapaya \
    Graimondii Tcacao \
    Osativa  OsativaKitaake \
    Bstacei Bdistachyon \
    ; do

    1>&2 echo "==> ${name}"

    pushd ~/data/alignment/JGI/${name} > /dev/null

    find ~/data/PhytozomeV12/${name}/annotation/ -name "*.gff3.gz" |
        grep "gene_exons" |
        xargs gzip -dcf > chr.gff
    
    # create anno.yml
    spanr gff --tag CDS chr.gff -o cds.yml
    spanr gff *.rm.gff -o repeat.yml
    spanr merge repeat.yml cds.yml -o anno.yml

    rm -f repeat.yml cds.yml

    popd  > /dev/null

done

```

## Athaliana

```bash
mkdir -p ~/data/alignment/JGI
cd ~/data/alignment/JGI

bash ~/data/alignment/JGI/jgi_fasta.sh Athaliana
```

Check assemblies with Ensembl version. Set `-l 100000` to avoid reverse matches.

**Assemblies matched.**

```bash
for chr in 1 2 3 4 5 ; do
    sparsemem -maxmatch -F -l 100000 -b -n -k 3 -threads 3 \
        ~/data/alignment/JGI/Athaliana/Chr${chr}.fa \
        ~/data/alignment/Ensembl/Atha/${chr}.fa
done

```

