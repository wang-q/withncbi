# HMM related resources

[TOC levels=1-3]: # " "
- [HMM related resources](#hmm-related-resources)
- [PFAM-A](#pfam-a)
- [TIGRFAM](#tigrfam)
- [40 single-copy genes](#40-single-copy-genes)


# PFAM-A

```bash
mkdir -p ~/data/HMM/PFAM

cd ~/data/HMM/PFAM

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
    wget -N --content-disposition ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/${basename}.gz
done

for basename in Pfam-A.hmm Pfam-A.hmm.dat active_site.dat; do
    echo "==> ${basename}"
    gzip -dcf ${basename}.gz > ${basename}
done

```

# TIGRFAM

```bash
mkdir -p ~/data/HMM/TIGRFAM

cd ~/data/HMM/TIGRFAM

wget -N --content-disposition ftp://ftp.jcvi.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz

mkdir -p HMM
tar xvfz TIGRFAMs_14.0_HMM.tar.gz --directory HMM

```

# 40 single-copy genes

Ref.:

1. Wu, D., Jospin, G. & Eisen, J. A. Systematic Identification of Gene Families for Use as “Markers”
   for Phylogenetic and Phylogeny-Driven Ecological Studies of Bacteria and Archaea and Their Major
   Subgroups. PLoS ONE 8, e77033 (2013).

2. Skennerton, C. T. et al. Phylogenomic analysis of Candidatus ‘Izimaplasma’ species: free-living
   representatives from a Tenericutes clade found in methane seeps. ISME J. 10, 2679–2692 (2016).

3. https://doi.org/10.6084/m9.figshare.722713.v1

4. `bacteria_and_archaea.tgz`: https://ndownloader.figshare.com/files/3093482

```bash
mkdir -p ~/data/HMM/40sg
cd ~/data/HMM/40sg

wget -N --content-disposition https://ndownloader.figshare.com/files/3093482

tar xvfz bacteria_and_archaea.tgz

```

