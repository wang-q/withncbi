# HMM related resources

- [HMM related resources](#hmm-related-resources)
    * [PFAM-A](#pfam-a)
    * [TIGRFAM](#tigrfam)
    * [40 single-copy genes](#40-single-copy-genes)

## PFAM-A

```shell script
mkdir -p ~/data/HMM/PFAM

cd ~/data/HMM/PFAM

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    proxychains4 wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
done

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done

```

## TIGRFAM

```shell script
mkdir -p ~/data/HMM/TIGRFAM

cd ~/data/HMM/TIGRFAM

proxychains4 wget -N --content-disposition ftp://ftp.jcvi.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz

mkdir -p HMM
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR02013.HMM
tar --directory HMM -xzvf TIGRFAMs_14.0_HMM.tar.gz TIGR00485.HMM

```

## 40 single-copy genes

Ref.:

1. Wu, D., Jospin, G. & Eisen, J. A. Systematic Identification of Gene Families for Use as “Markers”
   for Phylogenetic and Phylogeny-Driven Ecological Studies of Bacteria and Archaea and Their Major
   Subgroups. PLoS ONE 8, e77033 (2013).

2. Skennerton, C. T. et al. Phylogenomic analysis of Candidatus ‘Izimaplasma’ species: free-living
   representatives from a Tenericutes clade found in methane seeps. ISME J. 10, 2679–2692 (2016).

3. https://doi.org/10.6084/m9.figshare.722713.v1

4. `bacteria_and_archaea.tgz`: https://ndownloader.figshare.com/files/3093482

```shell script
mkdir -p ~/data/HMM/40scg
cd ~/data/HMM/40scg

curl -LO https://ndownloader.figshare.com/files/3093482

tar xvfz 3093482

```

