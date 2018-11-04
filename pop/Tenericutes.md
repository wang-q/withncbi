# Aligning various genera from Tenericutes


[TOC levels=1-3]: # " "
- [Aligning various genera from Tenericutes](#aligning-various-genera-from-tenericutes)
- [Phylum Tenericutes](#phylum-tenericutes)
- [Trichoderma: assembly](#trichoderma-assembly)
- [Count strains](#count-strains)
- [Find all RNase R](#find-all-rnase-r)
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
Paus_str_NZSb11
PAyel_AYWB
PBnap
Pbra_str_JR1
PCcor
Pcyel
PEpur
PIclo_str_MA1
PMyel_str_MW1
PNJer
POyel_OY_M
PRora
Psp_Vc33
PVwit_str_VAC
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

wc -l taxon/*

find taxon -maxdepth 1 -type f -not -name "*.replace.tsv" |
    xargs -i basename {} \
    > genus.list

```

| Order             | Genus           | Comments | Species | Strains |
|:------------------|:----------------|:---------|--------:|--------:|
| Acholeplasmatales | Acholeplasma    | 无胆甾原体 |      11 |      13 |
|                   | Phytoplasma     | 植原体    |         |      18 |
| Anaeroplasmatales | Anaeroplasma    |          |       1 |       1 |
|                   | Asteroleplasma  |          |       0 |       0 |
| Entomoplasmatales | Entomoplasma    |          |       6 |      10 |
|                   | Mesoplasma      |          |      11 |      15 |
|                   | Spiroplasma     | 螺原体    |      25 |      29 |
| Mycoplasmatales   | Mycoplasma      | 支原体    |      81 |     192 |
|                   | Ureaplasma      | 脲原体    |       4 |      25 |
| Unclassified      | Izimaplasma     | 独立生活   |       1 |       1 |
| Haloplasmatales   | Haloplasma      |          |       1 |       1 |
|                   | Inordinaticella |          |       1 |       1 |


# Find all RNase R

```bash
cd ~/data/alignment/Tenericutes

mkdir -p RnaseR

# 306
find ASSEMBLY -maxdepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    wc -l

# 303
find ASSEMBLY -type f -name "*_protein.faa.gz" |
    wc -l

find ASSEMBLY -type f -name "*_protein.faa.gz" |
    xargs gzip -dcf \
    > RnaseR/all.pro.fa

#find ASSEMBLY -type f -name "*_translated_cds.faa.gz" |
#    xargs gzip -dcf \
#    > RnaseR/all.tcds.fa
#
#cat RnaseR/all.tcds.fa |
#    grep "ribonuclease R"

# 280; deduped 207
faops some RnaseR/all.pro.fa \
    <(cat RnaseR/all.pro.fa |
        grep "ribonuclease R" |
        cut -d" " -f 1 |
        sed "s/^>//" |
        sort | uniq) \
    stdout |
    faops filter -u stdin stdout \
    > RnaseR/RnaseR.pro.fa

find ASSEMBLY -maxdepth 1 -type d |
    sort |
    grep 'ASSEMBLY/' |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        gzip -dcf {}/*_protein.faa.gz |
            grep "ribonuclease R" |
            (echo {} && cat)
        echo
    ' \
    > RnaseR/strains.txt

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
                    perl -nle '\''
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

# 280
#cat taxon/*.replace.tsv | wc -l

# extract sequences for each genus
for GENUS in $(cat genus.list); do
    echo "==> ${GENUS}"
    
    mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

    # avoid duplicated fasta headers
    faops some RnaseR/all.pro.fa taxon/${GENUS}.replace.tsv stdout |
        faops filter -u stdin ${mytmpdir}/${GENUS}.fa
    
    # avoid duplicated original names
    cat taxon/${GENUS}.replace.tsv |
        parallel --no-run-if-empty --linebuffer -k -j 1 "
            faops replace -s ${mytmpdir}/${GENUS}.fa <(echo {}) stdout
        " \
        > RnaseR/${GENUS}.pro.fa
        
    rm -fr ${mytmpdir}
done

# aligning with muscle
cat genus.list |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo "==> {}"
        
        muscle -quiet -in RnaseR/{}.pro.fa -out RnaseR/{}.aln.fa
    '

# newick trees
cat genus.list |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo "==> {}"
        
        FastTree -quiet RnaseR/{}.aln.fa > RnaseR/{}.aln.newick
    '

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

