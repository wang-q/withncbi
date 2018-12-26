# Aligning various genera from Tenericutes


[TOC levels=1-3]: # " "
- [Aligning various genera from Tenericutes](#aligning-various-genera-from-tenericutes)
- [Phylum Tenericutes](#phylum-tenericutes)
- [Trichoderma: assembly](#trichoderma-assembly)
    - [NCBI taxonomy](#ncbi-taxonomy)
- [Count strains](#count-strains)
- [Collect proteins](#collect-proteins)
- [Phylogenetics with 40 single-copy genes, *RpoB*, *EF-tu* and RNase_R](#phylogenetics-with-40-single-copy-genes-rpob-ef-tu-and-rnase_r)
    - [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    - [Create valid marker gene list](#create-valid-marker-gene-list)
    - [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
    - [Tweak the concat tree](#tweak-the-concat-tree)
    - [TIGR](#tigr)
- [RNase R](#rnase-r)
    - [Domain organisation](#domain-organisation)
        - [Human RRP44](#human-rrp44)
        - [Fly Dis3](#fly-dis3)
        - [Yeast RRP44](#yeast-rrp44)
        - [E. coli RNase R](#e-coli-rnase-r)
        - [E. coli RNase II](#e-coli-rnase-ii)
    - [RNase R domains](#rnase-r-domains)
    - [Scan every domains](#scan-every-domains)
    - [Stats of annotations and HMM models](#stats-of-annotations-and-hmm-models)
    - [Find all RNase R](#find-all-rnase-r)
    - [Tweak the tree of RNaseR](#tweak-the-tree-of-rnaser)
    - [RNB domain](#rnb-domain)
    - [Importin_rep domain](#importin_rep-domain)
- [Tenericutes: run](#tenericutes-run)


# Phylum Tenericutes

无壁菌门, 或称柔膜菌门, 无细胞壁

Ref.:

1. Skennerton, C. T. et al. Phylogenomic analysis of Candidatus ‘Izimaplasma’ species: free-living
   representatives from a Tenericutes clade found in methane seeps. ISME J. 10, 2679–2692 (2016).

2. Davis, J. J., Xia, F., Overbeek, R. A. & Olsen, G. J. Genomes of the class Erysipelotrichia
   clarify the firmicute origin of the class Mollicutes. Int. J. Syst. Evol. Microbiol. 63,
   2727–2741 (2013).

Key genera:

* Acholeplasmatales
    * *Acholeplasma*: 2147
    * *Candidatus Phytoplasma*: 33926

* Anaeroplasmatales
    * *Anaeroplasma*: 2086
    * *Asteroleplasma*: 2152

* Entomoplasmatales
    * *Entomoplasma*: 46238
    * *Mesoplasma*: 46239
    * *Spiroplasma*: 2132

* Mycoplasmatales
    * *Mycoplasma*: 2093
    * *Ureaplasma*: 2129
    * *Candidatus Hepatoplasma*: 295595

* *Candidatus Izimaplasma*: 1912503

* Haloplasmatales
    * *Haloplasma*: 471824
    * *Inordinaticella*: 1979191

Outgroup:

* Firmicutes
    * *Bacillus subtilis* subsp. subtilis str. 168: 224308
    * *Turicibacter sanguinis* PC909: 702450
    * Eubacterium limosum KIST612: 903814
    * Holdemania filiformis DSM 12042: 545696
    * Bulleidia extructa W1219: 679192
    * Solobacterium moorei F0204: 706433
    * Erysipelothrix rhusiopathiae str. Fujisawa: 650150
    * Erysipelothrix larvae: 1514105
    * Catenibacterium mitsuokai DSM 15897: 451640
    * Coprobacillus cateniformis: 100884
    * Clostridium acetobutylicum ATCC 824: 272562
    * Clostridium tetani E88: 212717
    * Clostridium botulinum A str. ATCC 3502: 413999

* Actinobacteria
    * Amycolatopsis mediterranei U32: 749927
    * Bifidobacterium adolescentis ATCC 15703: 367928
    * Corynebacterium glutamicum ATCC 13032: 196627
    * Mycobacterium tuberculosis H37Rv: 83332

*  Gammaproteobacteria
    * Escherichia coli str. K-12 substr. MG1655: 511145
    * Salmonella enterica subsp. enterica serovar Typhimurium str. LT2: 99287

Check NCBI pages:

* http://www.ncbi.nlm.nih.gov/assembly/?term=txid2093%5BOrganism:exp%5D
* http://www.ncbi.nlm.nih.gov/genome/?term=txid2093%5BOrganism:exp%5D

# Trichoderma: assembly

```bash
export RANK_NAME=Tenericutes

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (2147, 33926, 2086, 2152, 46238, 46239, 2132, 2093, 2129, 295595, 1912503, 471824, 1979191)
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND taxonomy_id in (224308, 702450, 903814, 545696, 679192, 706433, 650150, 1514105, 451640, 100884, 272562, 212717, 413999)
    " \
    >> raw.tsv

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND taxonomy_id in (749927, 367928, 196627, 83332, 511145, 99287)
    " \
    >> raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (2147, 33926, 2086, 2152, 46238, 46239, 2132, 2093, 2129, 295595, 1912503, 471824, 1979191)
    " \
    >> raw.tsv

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
    > ${RANK_NAME}.assembly.tsv

# comment out unneeded assembly levels

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

unset RANK_NAME

```

```bash
cd ~/data/alignment/Tenericutes

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/Tenericutes.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/Tenericutes.assembly.rsync.sh

bash ASSEMBLY/Tenericutes.assembly.collect.sh

```

## NCBI taxonomy

```bash
cd ~/data/alignment/Tenericutes

bp_taxonomy2tree.pl -e \
    -s "Acholeplasma" \
    -s "Anaeroplasma" \
    -s "Asteroleplasma" \
    -s "Entomoplasma" \
    -s "Haloplasma" \
    -s "Hepatoplasma" \
    -s "Inordinaticella" \
    -s "Izimaplasma" \
    -s "Mesoplasma" \
    -s "Mycoplasma" \
    -s "Phytoplasma" \
    -s "Spiroplasma" \
    -s "Ureaplasma" \
    -s "Bacillus subtilis" \
    -s "Bulleidia extructa" \
    -s "Catenibacterium mitsuokai" \
    -s "Coprobacillus cateniformis" \
    -s "Clostridium acetobutylicum" \
    -s "Clostridium botulinum" \
    -s "Clostridium tetani" \
    -s "Erysipelothrix larvae" \
    -s "Erysipelothrix rhusiopathiae" \
    -s "Eubacterium limosum" \
    -s "Holdemania filiformis" \
    -s "Solobacterium moorei" \
    -s "Turicibacter sanguinis" \
    -s "Amycolatopsis mediterranei" \
    -s "Bifidobacterium adolescentis" \
    -s "Corynebacterium glutamicum" \
    -s "Mycobacterium tuberculosis" \
    -s "Escherichia coli" \
    -s "Salmonella enterica" \
    > Tenericutes.newick

nw_display -w 600 -s Tenericutes.newick |
    rsvg-convert -o ~/Scripts/withncbi/pop/Tenericutes.png

```

![Tenericutes.png](../image/Tenericutes.png)

# Count strains

```bash
cd ~/data/alignment/Tenericutes

parallel --no-run-if-empty --linebuffer -k -j 4 '
    n_species=$(cat ASSEMBLY/Tenericutes.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        uniq |
        wc -l)
    
    n_strains=$(cat ASSEMBLY/Tenericutes.assembly.collect.csv |
        cut -d"," -f 2 |
        grep -v "Candidatus" |
        grep "{}" |
        cut -d" " -f 1,2 |
        sort |
        wc -l)
    
    printf "%s\t%d\t%d\n" {} ${n_species} ${n_strains}
    ' ::: Acholeplasma Entomoplasma Mesoplasma Spiroplasma Mycoplasma Ureaplasma

#Acholeplasma    11      12
#Entomoplasma    6       10
#Mesoplasma      11      15
#Spiroplasma     25      30
#Mycoplasma      86      197
#Ureaplasma      4       25

mkdir -p taxon

echo "An_bac" > taxon/Anaeroplasma
echo "Ha_cont_SSD_17B" > taxon/Haloplasma
echo "I_for" > taxon/Inordinaticella
echo "CH_crin_Av" > taxon/Hepatoplasma

cat <<EOF > taxon/Izimaplasma
CI_sp_HR1
CI_sp_HR2
CI_sp_Z
EOF

cat <<EOF > taxon/Phytoplasma
CP_Ast_AYWB
CP_Bra
CP_aura
CP_aus
CP_aus_Strawberry_NZSb11
CP_mali
CP_ory
CP_phoenici
CP_pru
CP_sol
CP_ziz
CP_Chrysanthemum_c
CP_Chrysanthemum_y
CP_Ech
CP_Ita_MA1
CP_Mai
CP_Mil_MW1
CP_Vac_VAC
CP_Whe
CP_New
CP_Oni_OY_M
CP_Pea_NTU2011
CP_Phy
CP_Poi_JR1
CP_Ric
EOF

parallel --no-run-if-empty --linebuffer -k -j 4 '
    cat ASSEMBLY/Tenericutes.assembly.collect.csv |
        cut -d"," -f 1,2 |
        grep "{}" |
        cut -d"," -f 1 \
        > taxon/{}
    ' ::: Acholeplasma Entomoplasma Mesoplasma Spiroplasma Mycoplasma Ureaplasma

# Misplaced in taxonomy tree
#echo "Ac_sp_CAG_878" >> taxon/Acholeplasma

#cat <<EOF >> taxon/Mycoplasma
#Mycop_sp_CAG_472
#Mycop_sp_CAG_611
#Mycop_sp_CAG_611_25_7
#Mycop_sp_CAG_776
#Mycop_sp_CAG_877
#Mycop_sp_CAG_956
#EOF

cat <<EOF > taxon/Actinobacteria
Am_med_U32
Bi_ado_ATCC_15703
Cor_glu_ATCC_13032
Mycob_tub_H37Rv
EOF

cat <<EOF > taxon/Gammaproteobacteria
Es_coli_K_12_MG1655
Sa_ente_Typhimurium_LT2
EOF

cat <<EOF > taxon/Clostridiales
Cl_ace_ATCC_824
Cl_bot_A_ATCC_3502
Cl_tet_E88
Eu_lim_KIST612
EOF

cat <<EOF > taxon/Erysipelotrichaceae
Bu_ext_W1219
Ca_mit_DSM_15897
Cop_cat
Er_lar
Er_rhu_Fujisawa
Ho_fil_DSM_12042
So_moo_F0204
EOF

cat <<EOF > taxon/Others
Ba_subt_subtilis_168
T_san_PC909
EOF

wc -l taxon/*

find taxon -maxdepth 1 -type f -not -name "*.replace.tsv" |
    xargs -i basename {} \
    > genus.list

# Omit strains without protein annotations
#Sp_Chol
#CP_Phy
#Mycop_moa_ATCC_27625
#Mycop_sp_Bg1
#Mycop_sp_Bg2
#Mycop_sp_U
for GENUS in $(cat genus.list); do
    for STRAIN in $(cat taxon/${GENUS}); do
        if ! compgen -G "ASSEMBLY/${STRAIN}/*_protein.faa.gz" > /dev/null; then
            echo ${STRAIN}
        fi
    done 
done \
    > omit.list

for GENUS in $(cat genus.list); do
    perl -nl -i -MPath::Tiny -e '
        BEGIN { 
            our %omit = map { $_ => 1 }
                        grep {/\S/}
                        path(q{omit.list})->lines({ chomp => 1});
        }
        print $_ if ! exists $omit{$_};
    ' taxon/${GENUS}
done

```

| Order             | Genus           | Comments           | Species | Strains |
|:------------------|:----------------|:-------------------|--------:|--------:|
| Acholeplasmatales |                 |                    |         |         |
|                   | Acholeplasma    | 无胆甾原体           |      11 |      13 |
|                   | Phytoplasma     | 植原体              |         |      25 |
| Anaeroplasmatales |                 |                    |         |         |
|                   | Anaeroplasma    |                    |       1 |       1 |
|                   | Asteroleplasma  |                    |         |         |
| Entomoplasmatales |                 | 虫原体              |         |         |
|                   | Entomoplasma    |                    |       6 |      10 |
|                   | Mesoplasma      |                    |      11 |      15 |
|                   | Spiroplasma     | 螺原体, 感染昆虫与植物 |      25 |      31 |
| Mycoplasmatales   |                 |                    |         |         |
|                   | Mycoplasma      | 支原体              |      86 |     207 |
|                   | Ureaplasma      | 脲原体              |       4 |      25 |
|                   | Hepatoplasma    |                    |       1 |       1 |
| Unclassified      |                 |                    |         |         |
|                   | Izimaplasma     | 独立生活             |         |       3 |
| Haloplasmatales   |                 |                    |         |         |
|                   | Haloplasma      |                    |       1 |       1 |
|                   | Inordinaticella |                    |       1 |       1 |

# Collect proteins

```bash
cd ~/data/alignment/Tenericutes

mkdir -p PROTEINS

# 352
find ASSEMBLY -maxdepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l

# 346
find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    xargs gzip -dcf \
    > PROTEINS/all.pro.fa

# ribonuclease
cat PROTEINS/all.pro.fa |
    grep "ribonuclease" |
    perl -nl -e 's/^>\w+\.\d+\s+//g; print' |
    perl -nl -e 's/\s+\[.+?\]$//g; print' |
    perl -nl -e 's/MULTISPECIES: //g; print' |
    sort |
    uniq -c -d |
    sort -nr

# methyltransferase
cat PROTEINS/all.pro.fa |
    grep "methyltransferase" |
    perl -nl -e 's/^>\w+\.\d+\s+//g; print' |
    perl -nl -e 's/\s+\[.+?\]$//g; print' |
    perl -nl -e 's/MULTISPECIES: //g; print' |
    sort |
    uniq -c -d |
    sort -nr

```

# Phylogenetics with 40 single-copy genes, *RpoB*, *EF-tu* and RNase_R

##  Find corresponding proteins by `hmmsearch`

* The `E_VALUE` was manually adjusted to 1e-20

* RpoB: TIGR02013
* EF-tu: TIGR00485

```bash
E_VALUE=1e-20

cd ~/data/alignment/Tenericutes

# example
#gzip -dcf ASSEMBLY/Aaxa_ATCC_25176/*_protein.faa.gz |
#    hmmsearch -E 1e-20 --domE 1e-20 --noali --notextw ~/data/HMM/40sg/bacteria_and_archaea_dir/BA00001.hmm - |
#    grep '>>' |
#    perl -nl -e '/>>\s+(\S+)/ and print $1'

# Find all genes
for marker in BA000{01..40}; do
    echo "==> marker [${marker}]"
    
    mkdir -p PROTEINS/${marker}
    
    for GENUS in $(cat genus.list); do
        echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/40scg/bacteria_and_archaea_dir/${marker}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ and printf qq{%s\t%s\n}, $1, $ENV{STRAIN};
                '
        done \
            > PROTEINS/${marker}/${GENUS}.replace.tsv
    done
    
    echo
done

for marker in TIGR02013 TIGR00485; do
    echo "==> marker [${marker}]"
    
    mkdir -p PROTEINS/${marker}
    
    for GENUS in $(cat genus.list); do
        echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw ~/data/HMM/TIGRFAM/HMM/${marker}.HMM - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '
                    />>\s+(\S+)/ or next;
                    $n = $1;
                    $s = $n;
                    $s =~ s/\.\d+//;
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                '
        done \
            > PROTEINS/${marker}/${GENUS}.replace.tsv
    done
    
    echo
done

```

## Create valid marker gene list

* `hmmsearch` may identify more than one copy for some marker genes.

    * BA00005: translation initiation factor IF-2
    * BA00008: signal recognition particle protein
    * BA00013: phenylalanine--tRNA ligase subunit beta

* Misidentified marker genes

    * BA00004: translation initiation factor EF-2

* Missing from many strains

    * BA00022: ribosomal protein L25/L23
    * BA00027: ribosomal protein L29
    * BA00032: tRNA pseudouridine synthase B
    * BA00035: Porphobilinogen deaminase
    * BA00038: phosphoribosylformylglycinamidine cyclo-ligase
    * BA00039: ribonuclease HII
    * BA00040: ribosomal protein L24

Compare proteins and strains.

```bash
cd ~/data/alignment/Tenericutes

for marker in BA000{01..03} BA000{06..07} BA000{09..12} BA000{14..21} BA000{23..26} BA000{28..31} BA000{33..34} BA000{36..37}; do
    echo ${marker}
done > marker.list

for marker in $(cat marker.list); do
    echo "==> marker [${marker}]"

    for GENUS in $(cat genus.list); do
        cat PROTEINS/${marker}/${GENUS}.replace.tsv |
            cut -f 2 |
            diff - taxon/${GENUS}
    done
    
    echo
done

```

## Align and concat marker genes to create species tree

```bash
cd ~/data/alignment/Tenericutes

# extract sequences 
for marker in $(cat marker.list) TIGR02013 TIGR00485; do
    echo "==> marker [${marker}]"

    for GENUS in $(cat genus.list); do
        echo "==> GENUS [${GENUS}]"

        mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

        # avoid duplicated fasta headers
        faops some PROTEINS/all.pro.fa PROTEINS/${marker}/${GENUS}.replace.tsv stdout |
            faops filter -u stdin ${mytmpdir}/${GENUS}.fa
        
        # avoid duplicated original names
        cat PROTEINS/${marker}/${GENUS}.replace.tsv |
            parallel --no-run-if-empty --linebuffer -k -j 1 "
                faops replace -s ${mytmpdir}/${GENUS}.fa <(echo {}) stdout
            " \
            > PROTEINS/${marker}/${GENUS}.pro.fa
            
        rm -fr ${mytmpdir}
    done
    
    echo
done

for marker in $(cat marker.list) TIGR02013 TIGR00485; do
    echo "==> marker [${marker}]"
    
    for GENUS in $(cat genus.list); do
        cat PROTEINS/${marker}/${GENUS}.pro.fa
    done \
        > PROTEINS/${marker}/${marker}.pro.fa
done

# aligning each markers with muscle
cat marker.list |
    parallel --no-run-if-empty --linebuffer -k -j 4 "
        echo '==> {}'
        
        muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    "

# concat marker genes
for marker in $(cat marker.list); do
    # sequences in one line
    faops filter -l 0 PROTEINS/${marker}/${marker}.aln.fa stdout
    
    # empty line for .fas
    echo
done \
    > PROTEINS/markers.aln.fas

# faspos names need full headers
#fasops names Phylo/markers.aln.fas -o stdout
cat PROTEINS/markers.aln.fas |
    grep "^>" |
    sed "s/^>//" |
    sort |
    uniq \
    > strains.list
fasops concat PROTEINS/markers.aln.fas strains.list -o PROTEINS/concat.aln.fa

FastTree -quiet PROTEINS/concat.aln.fa > PROTEINS/concat.aln.newick

```

## Tweak the concat tree

```bash
cd ~/data/alignment/Tenericutes

# reroot
nw_reroot PROTEINS/concat.aln.newick Am_med_U32 Es_coli_K_12_MG1655 > PROTEINS/concat.reroot.newick

# Check monophyly for genus
rm genus.monophyly.list genus.paraphyly.list genus.monophyly.map
for GENUS in $(cat genus.list); do
    NODE=$(
        nw_clade -m PROTEINS/concat.reroot.newick $(cat taxon/${GENUS}) |
            nw_stats -f l - |
            cut -f 3
    )
    
    if [[ "$NODE" ]]; then
        echo "${GENUS}" >> genus.monophyly.list
        cat taxon/${GENUS} |
            xargs -I{} echo "{} ${GENUS}___${NODE}" \
            >> genus.monophyly.map
    else
        echo "${GENUS}" >> genus.paraphyly.list
    fi
    
done

# Merge strains in genus to higher-rank
nw_rename PROTEINS/concat.reroot.newick genus.monophyly.map |
    nw_condense - \
    > PROTEINS/concat.map.newick

# strains in species
for GENUS in $(cat genus.paraphyly.list); do
    cat taxon/${GENUS}
done |
    perl -nl -e '/([[:alpha:]]+_[[:alpha:]]+)/ and print $1' |
    perl -nl -e '!/_sp$/ and print' |
    sort |
    uniq -d -c |
    perl -nla -e 'print qq{$F[1]\t$F[1]___$F[0]}' \
    > species.count.tsv

for GENUS in $(cat genus.paraphyly.list); do
    cat taxon/${GENUS}
done \
    > strains.paraphyly.list

# Check monophyly for species
rm species.monophyly.list species.paraphyly.list species.monophyly.map
cat species.count.tsv |
    perl -nl -MPath::Tiny -e '
        BEGIN {
            our @strains = 
                grep {/\S/}
                path(q{strains.paraphyly.list})->lines({ chomp => 1});
        }
        
        my @ns = split /\t/;
        my @sts = grep {/^$ns[0]/} @strains;
        
        my $cmd = q{nw_clade -m PROTEINS/concat.reroot.newick };
        $cmd .= " $_ " for @sts;
        $cmd .= " | nw_stats -f l - | cut -f 3";
        
        my $result = `$cmd`;
        if ($result) {
            print qq{$_ $ns[1]} for @sts;
            path(q{species.monophyly.list})->append($ns[0]);
        }
        else {
            path(q{species.paraphyly.list})->append($ns[0]);
        }
    ' \
    > species.monophyly.map

# Merge strains in species to higher-rank
nw_rename PROTEINS/concat.map.newick species.monophyly.map |
    nw_condense - \
    > PROTEINS/concat.map2.newick

```

## TIGR

```bash
cd ~/data/alignment/Tenericutes

parallel --no-run-if-empty --linebuffer -k -j 4 "
    echo '==> {}'
    
    muscle -quiet -in PROTEINS/{}/{}.pro.fa -out PROTEINS/{}/{}.aln.fa
    
    FastTree -quiet PROTEINS/{}/{}.aln.fa > PROTEINS/{}/{}.aln.newick
    " ::: TIGR02013 TIGR00485

```

# RNase R

## Domain organisation

### Human RRP44

* [Pfam](http://pfam.xfam.org/protein/RRP44_HUMAN)
* [UniProt](https://www.uniprot.org/uniprot/Q9Y2L1)
* [NCBI](https://www.ncbi.nlm.nih.gov/protein/NP_055768.3)
* Subcellular locationi: Nucleus

![RRP44_HUMAN](../image/RRP44_HUMAN.png)

### Fly Dis3

* [Pfam](http://pfam.xfam.org/protein/Q9VC93_DROME)
* [UniProt](https://www.uniprot.org/uniprot/Q9VC93)
* Subcellular locationi: Nucleus

![Dis3](../image/Q9VC93_DROME.png)

### Yeast RRP44

* [Pfam](http://pfam.xfam.org/protein/RRP44_YEAST)
* [UniProt](https://www.uniprot.org/uniprot/Q08162)
* [NCBI](https://www.ncbi.nlm.nih.gov/protein/NP_014621.1)
* Subcellular locationi: Nucleus, Mitochondrion

![RRP44_YEAST](../image/RRP44_YEAST.png)

### E. coli RNase R

* [Pfam](http://pfam.xfam.org/protein/RNR_ECOLI)
* [UniProt](https://www.uniprot.org/uniprot/P21499)
* [NCBI](https://www.ncbi.nlm.nih.gov/protein/NP_418600.4)

![RNR_ECOLI](../image/RNR_ECOLI.png)

### E. coli RNase II

* [Pfam](http://pfam.xfam.org/protein/RNB_ECOLI)
* [UniProt](https://www.uniprot.org/uniprot/P30850)
* [NCBI](https://www.ncbi.nlm.nih.gov/protein/NP_415802.1)

![RNB_ECOLI](../image/RNB_ECOLI.png)

## RNase R domains

* OB_RNB (PF08206)
* CSD2 (PF17876)
* RNB (PF00773)
* S1 (PF00575)

* HTH_12 (PF08461)
* RNase_II_C_S1 (PF18614)

* Importin_rep (PF18773)

* PIN_4 (PF13638)
* Rrp44_CSD1 (PF17216)
* OB_Dis3 (PF17849)
* Rrp44_S1 (PF17215)

* RNase_R: TIGR02063
* 3_prime_RNase: TIGR00358

```bash
cd ~/data/alignment/Tenericutes

mkdir -p DOMAINS/HMM
cd DOMAINS/HMM

for ID in PF08206 PF17876 PF00773 PF00575 \
    PF08461 PF18614 PF18773 \
    PF13638 PF17216 PF17849 PF17215; do
    wget -N --content-disposition http://pfam.xfam.org/family/${ID}/hmm
done

cp ~/data/HMM/TIGRFAM/HMM/TIGR02063.HMM TIGR02063.hmm
cp ~/data/HMM/TIGRFAM/HMM/TIGR00358.HMM TIGR00358.hmm

```

## Scan every domains

```bash
E_VALUE=1e-3

cd ~/data/alignment/Tenericutes

for domain in OB_RNB CSD2 RNB S1 \
    HTH_12 RNase_II_C_S1 Importin_rep \
    PIN_4 Rrp44_CSD1 OB_Dis3 Rrp44_S1 \
    TIGR02063 TIGR00358; do
    echo 1>&2 "==> domain [${domain}]"
        
    for GENUS in $(cat genus.list); do
        echo 1>&2 "==> GENUS [${GENUS}]"

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
    
    echo 1>&2
done

wc -l DOMAINS/*.replace.tsv
#   384 DOMAINS/OB_RNB.replace.tsv
#   282 DOMAINS/CSD2.replace.tsv
#   313 DOMAINS/RNB.replace.tsv
#   841 DOMAINS/S1.replace.tsv
#   127 DOMAINS/HTH_12.replace.tsv
#    90 DOMAINS/RNase_II_C_S1.replace.tsv
#     6 DOMAINS/Importin_rep.replace.tsv
#    14 DOMAINS/PIN_4.replace.tsv
#     5 DOMAINS/Rrp44_CSD1.replace.tsv
#     7 DOMAINS/OB_Dis3.replace.tsv
#     1 DOMAINS/Rrp44_S1.replace.tsv
#   390 DOMAINS/TIGR00358.replace.tsv
#   484 DOMAINS/TIGR02063.replace.tsv

# All proteins appeared
find DOMAINS/ -name "*.replace.tsv" |
    sort |
    parallel -j 1 'cut -f 2 {}' |
    sort -u \
    > DOMAINS/domains.tsv
wc -l DOMAINS/domains.tsv
#1111 DOMAINS/domains.tsv

# Status of domains
for domain in OB_RNB CSD2 RNB S1 \
    HTH_12 RNase_II_C_S1 Importin_rep \
    PIN_4 Rrp44_CSD1 OB_Dis3 Rrp44_S1 \
    TIGR02063 TIGR00358; do
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
#1111 lines, 14 fields

# Add header line
echo -e '#name\tOB_RNB\tCSD2\tRNB\tS1\tHTH_12\tRNase_II_C_S1\tImportin_rep\tPIN_4\tRrp44_CSD1\tOB_Dis3\tRrp44_S1\tTIGR02063\tTIGR00358' | 
    cat - DOMAINS/domains.tsv > temp && mv temp DOMAINS/domains.tsv

```

## Stats of annotations and HMM models

| Item               | Count |
|:-------------------|------:|
| " ribonuclease R " |   307 |
| " RNB "            |     6 |
| " RNase II "       |     1 |
| " RNase R "        |     2 |
| deduped            |   239 |
| OB_RNB             |   194 |
| CSD2               |   198 |
| RNB                |   224 |
| S1                 |   207 |

```bash
cd ~/data/alignment/Tenericutes

mkdir -p RNaseR

faops some PROTEINS/all.pro.fa \
    <(cat PROTEINS/all.pro.fa |
        grep -e " ribonuclease R " -e " RNB " -e " RNase II " -e " RNase R " |
        cut -d" " -f 1 |
        sed "s/^>//" |
        sort | uniq) \
    stdout |
    faops filter -u stdin stdout \
    > RNaseR/RNaseR.all.fa

cat PROTEINS/all.pro.fa |
    grep " ribonuclease R " |
    wc -l
cat PROTEINS/all.pro.fa |
    grep " RNB " |
    wc -l
cat PROTEINS/all.pro.fa |
    grep " RNase II " |
    wc -l
cat PROTEINS/all.pro.fa |
    grep " RNase R " |
    wc -l
cat RNaseR/RNaseR.all.fa |
    grep "^>" |
    wc -l
for domain in OB_RNB CSD2 RNB S1; do
    cat RNaseR/RNaseR.all.fa |
        grep "^>" |
        sed "s/^>//" |
        grep -Fx -f <(cut -f 1 DOMAINS/${domain}.replace.tsv) |
        wc -l
done

# Strains and RNase R
find ASSEMBLY -maxdepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        gzip -dcf {}/*_protein.faa.gz |
            grep -e " ribonuclease R " -e " RNB " -e " RNase II " -e " RNase R " |
            (echo {} && cat)
        echo
    ' \
    > RNaseR/strain_anno.txt

```

## Find all RNase R

```bash
cd ~/data/alignment/Tenericutes

# Find all genes
for GENUS in $(cat genus.list); do
    echo "==> GENUS [${GENUS}]"

    for STRAIN in $(cat taxon/${GENUS}); do
        gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
            grep -e " ribonuclease R " -e " RNB " -e " RNase II " -e " RNase R "|
            cut -d" " -f 1 |
            sed "s/^>//" |
            STRAIN=${STRAIN} perl -nl -MPath::Tiny -e '
                BEGIN {
                    our %seen = map {(split /\t/)[0] => 1} 
                        grep {/\S/}
                        path(q{DOMAINS/RNB.replace.tsv})->lines({ chomp => 1});
                }
                
                $n = $_;
                $s = $n;
                $s =~ s/\.\d+//;
                if (exists $seen{$n}) {
                    printf qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                }
                else {
                    printf STDERR qq{%s\t%s_%s\n}, $n, $ENV{STRAIN}, $s;
                }
            '
    done \
        > RNaseR/${GENUS}.replace.tsv
done

# 301
cat RNaseR/*.replace.tsv | wc -l

# extract sequences for each genus
for GENUS in $(cat genus.list); do
    echo "==> ${GENUS}"
    
    mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

    # avoid duplicated fasta headers
    faops some PROTEINS/all.pro.fa RNaseR/${GENUS}.replace.tsv stdout |
        faops filter -u stdin ${mytmpdir}/${GENUS}.fa
    
    # avoid duplicated original names
    cat RNaseR/${GENUS}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 1 "
            faops replace -s ${mytmpdir}/${GENUS}.fa <(echo {}) stdout
        " \
        > RNaseR/${GENUS}.pro.fa
        
    rm -fr ${mytmpdir}
done

# aligning with muscle
cat genus.list |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo "==> {}"
        
        muscle -quiet -in RNaseR/{}.pro.fa -out RNaseR/{}.aln.fa
    '

# newick trees
cat genus.list |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo "==> {}"
        
        FastTree -quiet RNaseR/{}.aln.fa > RNaseR/{}.aln.newick
    '

for GENUS in $(cat genus.list); do
    cat RNaseR/${GENUS}.pro.fa
done \
    > RNaseR/RNaseR.pro.fa

muscle -quiet -in RNaseR/RNaseR.pro.fa -out RNaseR/RNaseR.aln.fa
FastTree -quiet RNaseR/RNaseR.aln.fa > RNaseR/RNaseR.aln.newick

```

Annotated as RNase R but lacking RNB domains:

* OAL10424.1 Mycop_Chaemob_OAL10424
* AFO52271.1 Mycop_Chaemol_Purdue_AFO52271
* AEG72387.1 Mycop_haemof_Ohio2_AEG72387
* WP_112665413.1 Mycop_wen_WP_112665413
* NP_462407.3 Sa_ente_Typhimurium_LT2_NP_462407

## Tweak the tree of RNaseR

```bash
cd ~/data/alignment/Tenericutes

# reroot
nw_reroot RNaseR/RNaseR.aln.newick Am_med_U32_YP_003763410 > RNaseR/RNaseR.reroot.newick

nw_labels -I RNaseR/RNaseR.aln.newick > RNaseR.list

# strains in species
cat RNaseR.list |
    perl -nl -e '/([[:alpha:]]+_[[:alpha:]]+)/ and print $1' |
    perl -nl -e '!/_sp$/ and print' |
    sort |
    uniq -d -c |
    perl -nla -e 'print qq{$F[1]\t$F[1]___$F[0]}' \
    > species.count.tsv

# Check monophyly for species
rm species.monophyly.list species.paraphyly.list species.monophyly.map
cat species.count.tsv |
    perl -nl -MPath::Tiny -e '
        BEGIN {
            our @lists = 
                grep {/\S/}
                path(q{RNaseR.list})->lines({ chomp => 1});
        }
        
        my @ns = split /\t/;
        my @sts = grep {/^$ns[0]/} @lists;
        
        my $cmd = q{nw_clade -m RNaseR/RNaseR.reroot.newick };
        $cmd .= " $_ " for @sts;
        $cmd .= " | nw_stats -f l - | cut -f 3";
        
        my $result = `$cmd`;
        if ($result) {
            print qq{$_ $ns[1]} for @sts;
            path(q{species.monophyly.list})->append($ns[0]);
        }
        else {
            path(q{species.paraphyly.list})->append($ns[0]);
        }
    ' \
    > species.monophyly.map

# Merge strains in species to higher-rank
nw_rename RNaseR/RNaseR.reroot.newick species.monophyly.map |
    nw_condense - \
    > RNaseR/RNaseR.map.newick

```

## RNB domain

```bash
E_VALUE=1e-3

cd ~/data/alignment/Tenericutes

hmmsearch \
    -E ${E_VALUE} --domE ${E_VALUE} \
    -A RNaseR/RNB.sto DOMAINS/HMM/RNB.hmm RNaseR/RNaseR.pro.fa

esl-reformat fasta RNaseR/RNB.sto > RNaseR/RNB.pro.fa

muscle -quiet -in RNaseR/RNB.pro.fa -out RNaseR/RNB.aln.fa
FastTree -quiet RNaseR/RNB.aln.fa > RNaseR/RNB.aln.newick

# reroot
nw_reroot RNaseR/RNB.aln.newick Am_med_U32_YP_003763410/46-352 > RNaseR/RNB.reroot.newick

nw_labels -I RNaseR/RNB.aln.newick > RNB.list

# strains in species
cat RNB.list |
    perl -nl -e '/([[:alpha:]]+_[[:alpha:]]+)/ and print $1' |
    perl -nl -e '!/_sp$/ and print' |
    sort |
    uniq -d -c |
    perl -nla -e 'print qq{$F[1]\t$F[1]___$F[0]}' \
    > species.count.tsv

# Check monophyly for species
rm species.monophyly.list species.paraphyly.list species.monophyly.map
cat species.count.tsv |
    perl -nl -MPath::Tiny -e '
        BEGIN {
            our @lists = 
                grep {/\S/}
                path(q{RNB.list})->lines({ chomp => 1});
        }
        
        my @ns = split /\t/;
        my @sts = grep {/^$ns[0]/} @lists;
        
        my $cmd = q{nw_clade -m RNaseR/RNB.reroot.newick };
        $cmd .= " $_ " for @sts;
        $cmd .= " | nw_stats -f l - | cut -f 3";
        
        my $result = `$cmd`;
        if ($result) {
            print qq{$_ $ns[1]} for @sts;
            path(q{species.monophyly.list})->append($ns[0]);
        }
        else {
            path(q{species.paraphyly.list})->append($ns[0]);
        }
    ' \
    > species.monophyly.map

# Merge strains in species to higher-rank
nw_rename RNaseR/RNB.reroot.newick species.monophyly.map |
    nw_condense - \
    > RNaseR/RNB.map.newick

```


## Importin_rep domain

```bash
E_VALUE=1e-3

cd ~/data/alignment/Tenericutes

hmmsearch \
    -E ${E_VALUE} --domE ${E_VALUE} \
    -A RNaseR/Importin_rep.sto DOMAINS/HMM/Importin_rep.hmm RNaseR/RNaseR.pro.fa

esl-reformat fasta RNaseR/Importin_rep.sto > RNaseR/Importin_rep.pro.fa

muscle -quiet -in RNaseR/Importin_rep.pro.fa -out RNaseR/Importin_rep.aln.fa
FastTree -quiet RNaseR/Importin_rep.aln.fa > RNaseR/Importin_rep.aln.newick

```

# Tenericutes: run

```bash
$(brew --prefix)/Cellar/$(brew list --versions repeatmasker | sed 's/ /\//')/libexec/util/queryRepeatDatabase.pl \
    -stat -species Bacteria
```

* Rsync to hpcc

```bash
rsync -avP \
    ~/data/alignment/Tenericutes/ \
    wangq@202.119.37.251:data/alignment/Tenericutes

# rsync -avP wangq@202.119.37.251:data/alignment/Tenericutes/ ~/data/alignment/Tenericutes

```

`--perseq` for RefSeq_category Reference Genome assemblies.

```bash
cd ~/data/alignment/Tenericutes

# prep
egaz template \
    ASSEMBLY \
    --prep -o GENOMES \
    --perseq Mpne_M129 \
    --perseq Mflo_L1 \
    --perseq Mmyc_subsp_mycoides_SC_str_PG1 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--parallel 24"

bsub -q mpi -n 24 -J "T-0_prep" "bash GENOMES/0_prep.sh"

ls -t output.* | head -n 1 | xargs tail -f | grep "==>"

# gff
for n in Mpne_M129 Mflo_L1 Mmyc_subsp_mycoides_SC_str_PG1; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"
    
    gzip -dcf ${FILE_GFF} > GENOMES/${n}/chr.gff
done

# multi
egaz template \
    GENOMES/Mpne_M129 \
    $(find GENOMES -maxdepth 1 -type d -path "*/????*" | grep -v "Mpne_M129") \
    --multi -o multi/ \
    --rawphylo --parallel 24 -v

bsub -q mpi -n 24 -J "T-1_pair" "bash multi/1_pair.sh"
bsub -w "ended(T-1_pair)" \
    -q mpi -n 24 -J "T-2_rawphylo" "bash multi/2_rawphylo.sh"
bsub  -w "ended(T-2_rawphylo)" \
    -q mpi -n 24 -J "T-3_multi" "bash multi/3_multi.sh"

# multi_Pfal
egaz template \
    GENOMES/Pfal_3D7 \
    $(find GENOMES -maxdepth 1 -type d -path "*/????*" | grep "Pfal_" | grep -v "Pfal_3D7") \
    GENOMES/Prei_SY57 \
    --multi -o multi/ \
    --multiname multi_Pfal --tree multi/Results/multi.nwk --outgroup Prei_SY57 \
    --parallel 24 -v

bsub -q mpi -n 24 -J "T-3_multi" "bash multi/3_multi.sh"

# self
egaz template \
    GENOMES/Pfal_3D7 GENOMES/Pber_ANKA GENOMES/Pcha_chabaudi \
    GENOMES/Pcyn_strain_B GENOMES/Pkno_strain_H \
    --self -o self/ \
    --circos --parallel 24 -v

bsub -q mpi -n 24 -J "T-1_self" "bash self/1_self.sh"
bsub -w "ended(T-1_self)" \
    -q mpi -n 24 -J "T-3_proc" "bash self/3_proc.sh"
bsub  -w "ended(T-3_proc)" \
    -q mpi -n 24 -J "T-4_circos" "bash self/4_circos.sh"

```

