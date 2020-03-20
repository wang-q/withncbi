# Classifying Plasmids

[TOC levels=1-3]: # ""

- [Classifying Plasmids](#classifying-plasmids)
  - [NCBI RefSeq](#ncbi-refseq)
  - [MinHash to get non-redundant plasmids](#minhash-to-get-non-redundant-plasmids)


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
    parallel -j 4 '
        echo >&2 "==> {}"
        faops some refseq.fa {} stdout |
            mash sketch -k 21 -s 1000 -i -p 6 - -o {}.msh
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 4 '
        echo >&2 "==> {}"
        mash dist -p 6 {}.msh refseq.plasmid.k21s1000.msh > {}.tsv
    '

find job -maxdepth 1 -type f -name "[0-9]??" | sort |
    parallel -j 16 '
        cat {}.tsv |
            tsv-filter --ff-str-ne 1:2 --le 3:0.01
    ' \
    > redundant.tsv

```

