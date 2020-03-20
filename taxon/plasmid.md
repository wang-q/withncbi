# Classifying Plasmids

[TOC levels=1-3]: # ""

- [Classifying Plasmids](#classifying-plasmids)
  - [NCBI RefSeq](#ncbi-refseq)
  - [MinHash to get non-redundant plasmids](#minhash-to-get-non-redundant-plasmids)
  - [Grouping by MinHash](#grouping-by-minhash)


## NCBI RefSeq

```bash
mkdir -p ~/data/plasmid
cd ~/data/plasmid

rsync -avP ftp.ncbi.nlm.nih.gov::refseq/release/plasmid/ RefSeq/

gzip -dcf RefSeq/*.genomic.gbff.gz > genomic.gbff
perl ~/Scripts/withncbi/taxon/gb_taxon_locus.pl genomic.gbff > refseq_id_seq.csv
rm genomic.gbff

gzip -dcf RefSeq/plasmid.1.1.genomic.fna.gz |
    grep "^>" |
    head -n 5
#>NC_006130.1 Streptococcus pyogenes 71-724 plasmid pDN571, complete sequence
#>NC_004464.2 Citrobacter freundii plasmid pCTX-M3, complete sequence
#>NC_006427.1 Enterococcus faecium plasmid pJB01, complete sequence
#>NC_001370.1 Lactobacillus plantarum plasmid pC30il, complete sequence
#>NC_002810.1 Streptococcus mutans LM7 plasmid pLM7, complete sequence

faops n50 -S -C RefSeq/*.genomic.fna.gz
#N50     222278
#S       2072550889
#C       22389

gzip -dcf RefSeq/*.genomic.fna.gz > RefSeq/plasmid.fa

```

## MinHash to get non-redundant plasmids 

```bash
mkdir ~/data/plasmid/nr
cd ~/data/plasmid/nr

faops size ../RefSeq/plasmid.fa > refseq.sizes

tsv-filter refseq.sizes --le 2:2000 | wc -l
#2473

faops some ../RefSeq/plasmid.fa <(tsv-filter refseq.sizes --gt 2:2000) refseq.fa

cat refseq.fa |
    mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.plasmid.k21s1000.msh

# split
mkdir -p job
faops size refseq.fa |
    cut -f 1 |
    split -l 1000 -a 3 -d - job/

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some refseq.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.plasmid.k21s1000.msh > {}.tsv
    '

# distance < 0.01
find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 16 '
        cat {}.tsv |
            tsv-filter --ff-str-ne 1:2 --le 3:0.01
    ' \
    > redundant.tsv

head -n 5 redundant.tsv
#NZ_CP034776.1   NC_005249.1     0.000730741     0       970/1000
#NZ_CP034416.1   NC_005249.1     0.00580821      0       794/1000
#NZ_LR745046.1   NC_005249.1     0.0010072       0       959/1000
#NZ_LR745043.1   NC_005249.1     0.000656154     0       973/1000
#NZ_CP033694.1   NC_006323.1     0.00766986      0       741/1000

cat redundant.tsv | wc -l 
# 129384

cat redundant.tsv |
    perl -nla -F"\t" -MGraph::Undirected -e '
        BEGIN {
            our $g = Graph::Undirected->new;
        }
        
        $g->add_edge($F[0], $F[1]);
    
        END {
            for my $cc ( $g->connected_components ) {
                print join qq{\t}, sort @{$cc};
            }        
        }
    ' \
    > connected_components.tsv

cat connected_components.tsv |
    perl -nla -F"\t" -e 'printf qq{%s\n}, $_ for @F' \
    > components.list

wc -l connected_components.tsv components.list
#  2073 connected_components.tsv
#  9800 components.list

faops some -i refseq.fa components.list stdout > refseq.nr.fa
faops some refseq.fa <(cut -f 1 connected_components.tsv) stdout >> refseq.nr.fa

rm -fr job

```

## Grouping by MinHash

```bash
mkdir ~/data/plasmid/grouping
cd ~/data/plasmid/grouping

cat ../nr/refseq.nr.fa |
    mash sketch -k 21 -s 1000 -i -p 8 - -o refseq.nr.k21s1000.msh

# split
mkdir -p job
faops size ../nr/refseq.nr.fa |
    cut -f 1 |
    split -l 1000 -a 3 -d - job/

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        faops some ../nr/refseq.nr.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.nr.k21s1000.msh > {}.tsv
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 1 '
        cat {}.tsv |
            tsv-select -f 1-3
    ' \
    > dist_full.tsv

# 40 min
#cat dist_full.tsv |
#    Rscript -e '
#        library(readr);
#        library(tidyr);
#        library(ape);
#        pair_dist <- read_tsv(file("stdin"), col_names=F); 
#        tmp <- pair_dist %>%
#            pivot_wider( names_from = X2, values_from = X3, values_fill = list(X3 = 1.0) )
#        tmp <- as.matrix(tmp)
#        mat <- tmp[,-1]
#        rownames(mat) <- tmp[,1]
#        
#        dist_mat <- as.dist(mat)
#        clusters <- hclust(dist_mat, method = "ward.D2")
#        tree <- as.phylo(clusters) 
#        write.tree(phy=tree, file="tree.nwk")
#        
#        group <- cutree(clusters, h=0.5) # k=3
#        groups <- as.data.frame(group)
#        groups$ids <- rownames(groups)
#        rownames(groups) <- NULL
#        groups <- groups[order(groups$group), ]
#        write_tsv(groups, "groups.tsv")
#    '

# distance < 0.1
cat dist_full.tsv |
    tsv-filter --ff-str-ne 1:2 --le 3:0.1 \
    > connected.tsv

head -n 5 connected.tsv
#NZ_CP044448.1   NC_006994.1     0.0920479       0       78/1000
#NZ_CP030113.1   NC_002487.1     0.0974869       0       69/1000
#NC_017210.1     NC_009034.1     0.0955937       0       72/1000
#NZ_CP029748.1   NC_006671.1     0.0733545       0       120/1000
#NZ_CP011064.1   NC_006671.1     0.0579033       0       174/1000

cat connected.tsv | wc -l 
#337156

mkdir -p group
cat connected.tsv |
    perl -nla -F"\t" -MGraph::Undirected -MPath::Tiny -e '
        BEGIN {
            our $g = Graph::Undirected->new;
        }
        
        $g->add_edge($F[0], $F[1]);
    
        END {
            my @rare;
            my $serial = 1;
            my @ccs = $g->connected_components;
            @ccs = map { $_->[0] }
                sort { $b->[1] <=> $a->[1] }
                map { [ $_, scalar( @{$_} ) ] } @ccs;
            for my $cc ( @ccs ) {
                my $count = scalar @{$cc};
                if ($count < 50) {
                    push @rare, @{$cc};
                } 
                else {
                    path(qq{group/$serial.lst})->spew(map {qq{$_\n}} @{$cc});
                    $serial++;
                }
            }
            path(qq{group/00.lst})->spew(map {qq{$_\n}} @rare);
            
            path(qq{grouped.lst})->spew(map {qq{$_\n}} $g->vertices);
        }
    '

# get non-grouped
# this will no be divided to subgroups
faops some -i ../nr/refseq.nr.fa grouped.lst stdout |
    faops size stdin |
    cut -f 1 \
    > group/lonely.lst

wc -l group/*
#  2527 group/00.lst
#  4974 group/1.lst
#   198 group/2.lst
#   148 group/3.lst
#    93 group/4.lst
#    74 group/5.lst
#    73 group/6.lst
#    63 group/7.lst
#    52 group/8.lst
#    51 group/9.lst
#  3946 group/lonely.lst
# 12199 total

find group -maxdepth 1 -type f -name "[0-9]*.lst" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        
        faops some ../nr/refseq.nr.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
            
        mash dist -p 6 {}.msh {}.msh > {}.tsv
    '

find group -maxdepth 1 -type f -name "[0-9]*.lst.tsv" | sort |
    parallel -j 4 --line-buffer '
        echo >&2 "==> {}"
        
        cat {} |
            tsv-select -f 1-3 |
            Rscript -e '\''
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
                write.tree(phy=tree, file="{.}.tree.nwk")
                
                group <- cutree(clusters, h=0.5) # k=3
                groups <- as.data.frame(group)
                groups$ids <- rownames(groups)
                rownames(groups) <- NULL
                groups <- groups[order(groups$group), ]
                write_tsv(groups, "{.}.groups.tsv")
            '\''
    '

# subgroup
mkdir -p subgroup
cp group/lonely.lst subgroup/

find group -name "*.groups.tsv" | sort |
    parallel -j 1 -k '
        cat {} | sed -e "1d" | xargs -I[] echo "{/.}_[]"
    ' |
    sed -e 's/.lst.groups_/_/' |
    perl -na -F"\t" -MPath::Tiny -e '
        path(qq{subgroup/$F[0].lst})->append(qq{$F[1]});
    '
    
# append ccs
cat ../nr/connected_components.tsv |
    parallel -j 1 --colsep "\t" '
        file=$(rg -F -l  "{1}" subgroup)
        echo {} | tr "[:blank:]" "\n" >> ${file}
    '

wc -l subgroup/* |
    sort -nr |
    head -n 100

wc -l subgroup/* |
    perl -pe 's/^\s+//' |
    tsv-filter -d" " --ge 1:50 |
    wc -l
#76

wc -l subgroup/* |
    perl -pe 's/^\s+//' |
    tsv-filter -d" " --le 1:2 |
    wc -l
#249

```

