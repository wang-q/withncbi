# Aligning various genera from Tenericutes


[TOC levels=1-3]: # " "
- [Aligning various genera from Tenericutes](#aligning-various-genera-from-tenericutes)
- [Phylum Tenericutes](#phylum-tenericutes)
- [Trichoderma: assembly](#trichoderma-assembly)
- [Count strains](#count-strains)
- [Find all RNase R](#find-all-rnase-r)
- [Phylogenetics with 40 single-copy genes](#phylogenetics-with-40-single-copy-genes)
    - [Find corresponding proteins by `hmmsearch`](#find-corresponding-proteins-by-hmmsearch)
    - [Create valid marker gene list](#create-valid-marker-gene-list)
    - [Align and concat marker genes to create species tree](#align-and-concat-marker-genes-to-create-species-tree)
- [Tenericutes: run](#tenericutes-run)


# Phylum Tenericutes

柔膜菌门, 无细胞壁

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

* *Candidatus Izimaplasma*: 1912503

* Haloplasmatales
    * *Haloplasma*: 471824
    * *Inordinaticella*: 1979191

Ref.:

1. Skennerton, C. T. et al. Phylogenomic analysis of Candidatus ‘Izimaplasma’ species: free-living
   representatives from a Tenericutes clade found in methane seeps. ISME J. 10, 2679–2692 (2016).

2. Davis, J. J., Xia, F., Overbeek, R. A. & Olsen, G. J. Genomes of the class Erysipelotrichia
   clarify the firmicute origin of the class Mollicutes. Int. J. Syst. Evol. Microbiol. 63,
   2727–2741 (2013).

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

Check NCBI pages

* http://www.ncbi.nlm.nih.gov/assembly/?term=txid2093%5BOrganism%3Aexp
* http://www.ncbi.nlm.nih.gov/genome/?term=txid2093%5BOrganism:exp%5D

# Trichoderma: assembly

```bash
export RANK_NAME=Tenericutes

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (2147, 33926, 2086, 2152, 46238, 46239, 2132, 2093, 2129, 1912503, 471824, 1979191)
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND taxonomy_id in (224308, 702450, 903814, 545696, 679192, 706433, 650150, 1514105, 451640, 100884, 272562, 212717, 413999)
    " \
    >> raw.tsv

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND taxonomy_id in (83332, 196627, 749927, 367928)
    " \
    >> raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT 
        organism_name, species, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (2147, 33926, 2086, 2152, 46238, 46239, 2132, 2093, 2129, 1912503, 471824, 1979191)
    " \
    >> raw.tsv

cat raw.tsv |
    (echo -e '#name\tftp_path\torganism\tassembly_level' && cat ) |
    perl -nl -a -F"\t" -e '
        BEGIN{my %seen}; 
        /^#/ and print and next;
        /^organism_name/i and next;
        $n = $F[0];
        $rx = quotemeta $F[1];
        $n =~ s/$rx\s*//;
        $n =~ s/\s+$//;
        $n =~ s/\W+/_/g;
        @O = split(/ /, $F[1]);
        $name = substr($O[0],0,1) . substr($O[1],0,3);
        $name .= q{_} . $n if $n;
        $name =~ s/\W+/_/g;
        $name =~ s/_+/_/g;
        $seen{$name}++;
        $seen{$name} > 1 and next;
        printf qq{%s\t%s\t%s\t%s\n}, $name, $F[2], $F[1], $F[3];
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
#Anaeroplasma    1       1
#Entomoplasma    6       10
#Mesoplasma      11      15
#Spiroplasma     25      28
#Mycoplasma      81      183
#Ureaplasma      4       25

mkdir -p taxon

echo "Abac" > taxon/Anaeroplasma
echo "CIsp_HR1" > taxon/Izimaplasma
echo "Hcon_SSD_17B" > taxon/Haloplasma
echo "Ifor" > taxon/Inordinaticella

cat <<EOF > taxon/Phytoplasma
Paus
Paus_NZSb11
PAyel_AYWB
PBnap
Pbra_JR1
PCcor
Pcyel
PEpur
PIclo_MA1
PMyel_MW1
PNJer
POyel_OY_M
PRora
Psp_Vc33
PVwit_VAC
PWblu
Pwit_NTU2011
PZmay
EOF

parallel --no-run-if-empty --linebuffer -k -j 4 '
    cat ASSEMBLY/Tenericutes.assembly.collect.csv |
        cut -d"," -f 1,2 |
        grep "{}" |
        cut -d"," -f 1 \
        > taxon/{}
    ' ::: Acholeplasma Entomoplasma Mesoplasma Spiroplasma Mycoplasma Ureaplasma

echo "Asp_878" >> taxon/Acholeplasma
cat <<EOF >> taxon/Mycoplasma
Msp_472
Msp_611
Msp_611_25_7
Msp_776
Msp_877
Msp_956
EOF

cat <<EOF > taxon/Outgroup
Bsub_subtilis_168
Bext_W1219
Cmit_DSM_15897
Cace_ATCC_824
Cbot_A_ATCC_3502
Ctet_E88
Ccat
Elar
Erhu_Fujisawa
Elim_KIST612
Hfil_DSM_12042
Smoo_F0204
Tsan_PC909
Amed_U32
Bado_ATCC_15703
Cglu_ATCC_13032
Mtub_H37Rv
EOF

wc -l taxon/*

find taxon -maxdepth 1 -type f -not -name "*.replace.tsv" |
    xargs -i basename {} \
    > genus.list

# Omit strains without protein annotations
#CShol
#Psp_Vc33
#Mmoa_ATCC_27625
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
| Acholeplasmatales | Acholeplasma    | 无胆甾原体           |      11 |      13 |
|                   | Phytoplasma     | 植原体              |         |      18 |
| Anaeroplasmatales | Anaeroplasma    |                    |       1 |       1 |
|                   | Asteroleplasma  |                    |       0 |       0 |
| Entomoplasmatales | Entomoplasma    |                    |       6 |      10 |
| 虫原体             | Mesoplasma      |                    |      11 |      15 |
|                   | Spiroplasma     | 螺原体, 感染昆虫与植物 |      25 |      29 |
| Mycoplasmatales   | Mycoplasma      | 支原体              |      81 |     192 |
|                   | Ureaplasma      | 脲原体              |       4 |      25 |
| Unclassified      | Izimaplasma     | 独立生活             |       1 |       1 |
| Haloplasmatales   | Haloplasma      |                    |       1 |       1 |
|                   | Inordinaticella |                    |       1 |       1 |


# Find all RNase R

```bash
cd ~/data/alignment/Tenericutes

mkdir -p RNaseR

# 319
find ASSEMBLY -maxdepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l

# 316
find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    xargs gzip -dcf \
    > RNaseR/all.pro.fa

#find ASSEMBLY -type f -name "*_translated_cds.faa.gz" |
#    xargs gzip -dcf \
#    > RNaseR/all.tcds.fa
#
#cat RNaseR/all.tcds.fa |
#    grep "ribonuclease R"

# 293; deduped 220
faops some RNaseR/all.pro.fa \
    <(cat RNaseR/all.pro.fa |
        grep "ribonuclease R" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        sort | uniq) \
    stdout |
    faops filter -u stdin stdout \
    > RNaseR/RNaseR.pro.fa

cat RNaseR/all.pro.fa |
    grep "ribonuclease R" |
    wc -l
cat RNaseR/RNaseR.pro.fa |
    grep "^>" |
    wc -l

muscle -quiet -in RNaseR/RNaseR.pro.fa -out RNaseR/RNaseR.aln.fa
FastTree -quiet RNaseR/RNaseR.aln.fa > RNaseR/RNaseR.aln.newick

find ASSEMBLY -maxdepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        gzip -dcf {}/*_protein.faa.gz |
            grep "ribonuclease R" |
            (echo {} && cat)
        echo
    ' \
    > RNaseR/strains.txt

for GENUS in $(cat genus.list); do
    cat taxon/${GENUS} |
        parallel --no-run-if-empty --linebuffer -k -j 4 '
            result=$(
                gzip -dcf ASSEMBLY/{}/*_protein.faa.gz |
                    grep "ribonuclease R" |
                    cut -d" " -f 1 |
                    sed "s/^>//"
            )
            if [ "$result" ]; then
                echo $result |
                    perl -nl -e '\''
                        @ns = split /\s+/;
                        for $n (@ns) {
                            $s = $n;
                            $s =~ s/\.\d+//;
                            printf qq{%s\t%s_%s\n}, $n, {}, $s;
                        }
                    '\''
            fi
        ' \
        > taxon/${GENUS}.replace.tsv
done

# 293
cat taxon/*.replace.tsv | wc -l

# extract sequences for each genus
for GENUS in $(cat genus.list); do
    echo "==> ${GENUS}"
    
    mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

    # avoid duplicated fasta headers
    faops some RNaseR/all.pro.fa taxon/${GENUS}.replace.tsv stdout |
        faops filter -u stdin ${mytmpdir}/${GENUS}.fa
    
    # avoid duplicated original names
    cat taxon/${GENUS}.replace.tsv |
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

```

# Phylogenetics with 40 single-copy genes

##  Find corresponding proteins by `hmmsearch`

The `E_VALUE` was manually adjusted to 1e-20.

```bash
E_VALUE=1e-20

cd ~/data/alignment/Tenericutes

# example
#gzip -dcf ASSEMBLY/Aaxa_ATCC_25176/*_protein.faa.gz |
#    hmmsearch -E 1e-20 --domE 1e-20 Phylo/bacteria_and_archaea_dir/BA00001.hmm - |
#    grep '>>' |
#    perl -nl -e '/>>\s+(\S+)/ and print $1'

# Find all marker genes
for marker in BA000{01..40}; do
    echo "==> marker [${marker}]"
    
    mkdir -p Phylo/${marker}
    
    for GENUS in $(cat genus.list); do
        echo "==> GENUS [${GENUS}]"

        for STRAIN in $(cat taxon/${GENUS}); do
            gzip -dcf ASSEMBLY/${STRAIN}/*_protein.faa.gz |
                hmmsearch -E ${E_VALUE} --domE ${E_VALUE} Phylo/bacteria_and_archaea_dir/${marker}.hmm - |
                grep '>>' |
                STRAIN=${STRAIN} perl -nl -e '/>>\s+(\S+)/ and printf qq{%s\t%s\n}, $1, $ENV{STRAIN}'
        done \
            > Phylo/${marker}/${GENUS}.replace.tsv
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
        cat Phylo/${marker}/${GENUS}.replace.tsv |
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
for marker in $(cat marker.list); do
    echo "==> marker [${marker}]"

    for GENUS in $(cat genus.list); do
        echo "==> GENUS [${GENUS}]"

        mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

        # avoid duplicated fasta headers
        faops some RNaseR/all.pro.fa Phylo/${marker}/${GENUS}.replace.tsv stdout |
            faops filter -u stdin ${mytmpdir}/${GENUS}.fa
        
        # avoid duplicated original names
        cat Phylo/${marker}/${GENUS}.replace.tsv |
            parallel --no-run-if-empty --linebuffer -k -j 1 "
                faops replace -s ${mytmpdir}/${GENUS}.fa <(echo {}) stdout
            " \
            > Phylo/${marker}/${GENUS}.pro.fa
            
        rm -fr ${mytmpdir}
    done
    
    echo
done

for marker in $(cat marker.list); do
    echo "==> marker [${marker}]"
    
    for GENUS in $(cat genus.list); do
        cat Phylo/${marker}/${GENUS}.pro.fa
    done \
        > Phylo/${marker}/${marker}.pro.fa
done

# aligning each markers with muscle
cat marker.list |
    parallel --no-run-if-empty --linebuffer -k -j 4 "
        echo '==> {}'
        
        muscle -quiet -in Phylo/{}/{}.pro.fa -out Phylo/{}/{}.aln.fa
    "

# concat marker genes
for marker in $(cat marker.list); do
    # sequences in one line
    faops filter -l 0 Phylo/${marker}/${marker}.aln.fa stdout
    
    # empty line for .fas
    echo
done \
    > Phylo/markers.aln.fas

# faspos names need full headers
#fasops names Phylo/markers.aln.fas -o stdout
cat Phylo/markers.aln.fas |
    grep "^>" |
    sed "s/^>//" |
    sort |
    uniq \
    > strains.list
fasops concat Phylo/markers.aln.fas strains.list -o Phylo/concat.aln.fa

FastTree -quiet Phylo/concat.aln.fa > Phylo/concat.aln.newick

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

