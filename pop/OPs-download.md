# Downloading steps for each groups

Less detailed than Trichoderma in [README.md](README.md), but include examples for genomes out of
WGS, which usually in better assembling levels.

[TOC]: # ""

- [*Arabidopsis* 19 genomes](#arabidopsis-19-genomes)
- [*Orazy sativa* Japonica 24 genomes](#orazy-sativa-japonica-24-genomes)
- [*Drosophila* Population Genomics Project (dpgp)](#drosophila-population-genomics-project-dpgp)
- [Primates](#primates)
- [Human individuals from Simons project](#human-individuals-from-simons-project)
- [*Caenorhabditis elegans*](#caenorhabditis-elegans)
- [Mouse](#mouse)
- [Currently not used](#currently-not-used)

## *Arabidopsis* 19 genomes

1. Sources.

   * [Project page](http://mus.well.ox.ac.uk/19genomes/)
   * [Download page](http://mus.well.ox.ac.uk/19genomes/fasta/)
   * [Paper](http://www.nature.com/nature/journal/v477/n7365/full/nature10414.html)

2. Download with my web page crawler.

   ```bash
   mkdir -p ~/data/alignment/others
   cd ~/data/alignment/others

   perl ~/Scripts/download/list.pl -u http://mus.well.ox.ac.uk/19genomes/fasta/
   perl ~/Scripts/download/download.pl -i 19genomes_fasta.yml -a
   aria2c -x 9 -s 3 -c -i /home/wangq/data/alignment/others/19genomes_fasta.yml.txt
   ```

3. *A. thaliana* and *A. lyrata* from ensembl genomes.

   [ensembl_82.yml](ensembl/ensembl_82.yml) creates `Atha` and `Alyr` in `~/data/alignment/Ensembl`.

## *Orazy sativa* Japonica 24 genomes

1. Sources.

   * [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP003189)
   * [Paper](http://www.nature.com/nbt/journal/v30/n1/full/nbt.2050.html)

   Mapping strategy in [here](https://github.com/wang-q/sra/blob/master/japonica24_seq.pl).

2. 23 Japonica rices restore from previous .2bit files.

   I've used bwa-gatk pipeline to generate 23 Japonica rice genomes.

   Reference assembly of nipponbare in that time was MSU6, now is IRGSP-1.0. But I don't want to
   waste time to rebuild and RepeatMasker all sequences.

   ```bash
   find ~/data/alignment/rice/ -name "*.2bit" \
       | grep -v "_65" \
       | parallel basename {//} \
       | sort

   mkdir -p ~/data/alignment/others/japonica24
   cd ~/data/alignment/others/japonica24

   for d in IRGC11010 IRGC1107 IRGC12793 IRGC17757 IRGC2540 IRGC26872 IRGC27630 IRGC31856 IRGC32399 IRGC328 IRGC38698 IRGC38994 IRGC418 IRGC43325 IRGC43675 IRGC50448 IRGC55471 IRGC66756 IRGC8191 IRGC8244 IRGC9060 IRGC9062 RA4952
   do
       twoBitToFa ~/data/alignment/rice/${d}/chr.2bit ${d}.fa;
   done
   ```

3. nipponbare and 9311 from ensembl genomes.

   [ensembl_82.yml](ensembl/ensembl_82.yml) creates `OsatJap` and `OsatInd` in
   `~/data/alignment/Ensembl`.

## *Drosophila* Population Genomics Project (dpgp)

1. Sources.

   * [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP005599)
   * [Paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003080)

   Mapping strategy in [here](https://github.com/wang-q/sra/blob/master/dpgp_seq.pl).

2. 21 genomes restore from previous .2bit files.

   ```bash
   find ~/data/alignment/dpgp/ -name "*.2bit" \
       | grep -v "_65" \
       | parallel basename {//} \
       | sort

   mkdir -p ~/data/alignment/others/dpgp
   cd ~/data/alignment/others/dpgp

   for d in CK1 CO15N CO2 CO4N ED10N EZ5N FR207 FR217 FR229 FR361 GA125 GA129 GA130 GA132 GA185 GU10 KN6 KR39 KT1 NG3N RC1 RG15 SP254 TZ8 UG7 UM526 ZI268 ZL130 ZO12 ZS37
   do
       twoBitToFa ~/data/alignment/dpgp/${d}/chr.2bit ${d}.fa;
   done
   ```

3. Dmel and Dsim from ensembl genomes.

   [ensembl_82.yml](ensembl/ensembl_82.yml) creates `Dmel` and `Dsim` in
   `~/data/alignment/Ensemble`.

## Primates

1. Guild tree

   ```bash
   cd ~/data/alignment
   wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/hg19.100way.commonNames.nh
   tree_doctor hg19.100way.commonNames.nh --newick  \
       --prune-all-but \
       Human,Chimp,Gorilla,Orangutan,Gibbon,Rhesus,Crab_eating_macaque,Baboon,Green_monkey,Marmoset,Squirrel_monkey,Bushbaby,Chinese_tree_shrew \
       > primates_13way.nwk

   ```

2. All from ensembl.

   [ensembl_82.yml](ensembl/ensembl_82.yml) creates `Human`, `Chimp`, `Gorilla`, `Orangutan` and
   `Rhesus` in `~/data/alignment/Ensembl`.

## Human individuals from Simons project

1. Sources

   * https://www.simonsfoundation.org/life-sciences/simons-genome-diversity-project-dataset/
   * Download with Globus
     * https://www.globus.org/xfer/StartTransfer?origin=simonsfoundation%23sf_data/SGVP/Whole_Genome/Public/HMS_Reich/&source_id=d20e60cc-6d04-11e5-ba46-22000b92c6ec&source_path=/SGVP/Whole_Genome/Public/HMS_Reich/
   * The data were mapped to hg19/GRCh37.
     * https://www.simonsfoundation.org/life-sciences/simons-genome-diversity-project/

2. The data are experimentally phased. Pick one from diploid genomes.

   ```bash
   mkdir -p ~/data/alignment/others/simons-phase0
   tree ~/data/alignment/others/simons-phase0
   ```

   ```text
   ~/data/alignment/others/simons-phase0
   ├── HGDP00456.phased.0.Q40.fa.gz
   ├── HGDP00521.phased.0.Q40.fa.gz
   ├── HGDP00542.phased.0.Q40.fa.gz
   ├── HGDP00665.phased.0.Q40.fa.gz
   ├── HGDP00778.phased.0.Q40.fa.gz
   ├── HGDP00927.phased.0.Q40.fa.gz
   ├── HGDP00998.phased.0.Q40.fa.gz
   ├── HGDP01029.phased.0.Q40.fa.gz
   ├── HGDP01284.phased.0.Q40.fa.gz
   ├── HGDP01307.phased.0.Q40.fa.gz
   └── MIXE.phased.0.Q40.fa.gz
   ```

## *Caenorhabditis elegans*

There are no suitable outgroups for *C. elegans*.

http://hgdownload.soe.ucsc.edu/goldenPath/ce10/multiz7way/ce10.commonNames.7way.nh

1. 40 wild strains from cele_mmp.

   Mapping strategy in [here](https://github.com/wang-q/sra/blob/master/cele_mmp_seq.pl).

   ```bash
   mkdir -p ~/data/alignment/others/cele
   cd ~/data/alignment/others/cele

   find ~/data/dna-seq/cele_mmp/ -name "*.vcf.fasta" \
       | parallel -j 1 cp {} .
   ```

2. Reference strain N2 from ensembl genomes

   [ensembl_82.yml](ensembl/ensembl_82.yml) creates `Cele` in `~/data/alignment/Ensembl`.

## Mouse

1. Sources.

   * [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP000927)
   * [Project website](http://www.sanger.ac.uk/science/data/mouse-genomes-project)
   * [Paper](http://dx.doi.org/10.1038/nature10413)

2. Download masked de novo assemblies from sanger-mouse ftp server.

   ```bash
   mkdir -p ~/data/alignment/others/sanger-mouse
   cd ~/data/alignment/others/sanger-mouse
   wget -m ftp://ftp-mouse.sanger.ac.uk/REL-1509-Assembly --accept "*.fa.masked.gz,README" .
   mv ftp-mouse.sanger.ac.uk/REL-1509-Assembly/* .
   rm -fr ftp-mouse.sanger.ac.uk

   # rsync -avP wangq@45.79.80.100:data/alignment/others/sanger-mouse/ ~/data/alignment/others/sanger-mouse
   ```

3. Reference strain C57BL/6J (GRCm38) and rat from ensembl.

   [ensembl_82.yml](ensembl/ensembl_82.yml) creates `Mouse` and `Rat` in `~/data/alignment/Ensembl`.

## Currently not used

* Genus *Drosophila*

  * Phylogenetic tree

    ```bash
    # download from http://hgdownload.soe.ucsc.edu/goldenPath/dm6/multiz27way/
    ~/share/phast/bin/tree_doctor dm6.27way.scientificNames.nh --newick \
        --rename "Drosophila_melanogaster -> Dmel ; Drosophila_simulans -> Dsim ; Drosophila_sechellia -> Dsech" \
        > temp1.nwk

    ~/share/phast/bin/tree_doctor temp1.nwk --newick \
        --rename "Drosophila_yakuba -> Dyak ; Drosophila_erecta -> Dere" \
        > temp2.nwk

    ~/share/phast/bin/tree_doctor temp2.nwk --newick \
        --rename "Drosophila_pseudoobscura_pseudoobscura -> Dpse ; Drosophila_persimilis -> Dper" \
        > temp3.nwk

    ~/share/phast/bin/tree_doctor temp3.nwk --newick \
        --prune-all-but Dmel,Dsim,Dsech,Dyak,Dere,Dpse,Dper \
        > fly_7way.nwk

    rm temp[0-9].nwk
    ```

  * taxon

    ```bash
    --download 'name=Dmel;taxon=7227;sciname=Drosophila melanogaster' \
    --download 'name=Dsim;taxon=7240;sciname=Drosophila simulans' \
    --download 'name=Dsech;taxon=7238;sciname=Drosophila sechellia' \
    --download 'name=Dyak;taxon=7245;sciname=Drosophila yakuba' \
    --download 'name=Dere;taxon=7220;sciname=Drosophila erecta' \
    --download 'name=Dpse;taxon=7237;sciname=Drosophila pseudoobscura' \
    --download 'name=Dper;taxon=7234;sciname=Drosophila persimilis simulans' \
    ```

* DGRP

  * taxon

    ```perl
    my @data = (
        { taxon => 900501, name => "DGRP-138", coverage => 34, },
        { taxon => 900502, name => "DGRP-176", coverage => 38.48, },
        { taxon => 900503, name => "DGRP-181", coverage => 30.66, },
        { taxon => 900504, name => "DGRP-208", coverage => 34.51, },
        { taxon => 900505, name => "DGRP-321", coverage => 41.75, },
        { taxon => 900506, name => "DGRP-332", coverage => 31.53, },
        { taxon => 900507, name => "DGRP-375", coverage => 41.91, },
        { taxon => 900508, name => "DGRP-38",  coverage => 34.1, },
        { taxon => 900509, name => "DGRP-380", coverage => 36.73, },
        { taxon => 900510, name => "DGRP-391", coverage => 47.62, },
        { taxon => 900511, name => "DGRP-40",  coverage => 41.27, },
        { taxon => 900512, name => "DGRP-406", coverage => 30.3, },
        { taxon => 900513, name => "DGRP-443", coverage => 33.35, },
        { taxon => 900514, name => "DGRP-517", coverage => 45.97, },
        { taxon => 900515, name => "DGRP-57",  coverage => 38.28, },
        { taxon => 900516, name => "DGRP-727", coverage => 33.64, },
        { taxon => 900517, name => "DGRP-738", coverage => 31.74, },
        { taxon => 900518, name => "DGRP-757", coverage => 32.87, },
        { taxon => 900519, name => "DGRP-852", coverage => 40.42, },
        { taxon => 900520, name => "DGRP-897", coverage => 32.81, },
    );

    ```

* Primates

  * taxon

    ```bash
    --download 'name=Gibbon;taxon=61853;sciname=Nomascus leucogenys;coverage=5.6x sanger' \
    --download 'name=Marmoset;taxon=9483;sciname=Callithrix jacchus;coverage=6x sanger' \
    --download 'name=Tarsier;taxon=9478;sciname=Tarsius syrichta;coverage=1.82x sanger' \
    --download 'name=Lemur;taxon=30608;sciname=Microcebus murinus;coverage=1.93x sanger' \
    --download 'name=Bushbaby;taxon=30611;sciname=Otolemur garnettii;coverage=2x sanger' \
    ```

