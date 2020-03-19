# Aligning various genera from Protists

Less detailed than *Trichoderma* in
[README.md](https://github.com/wang-q/withncbi/blob/master/pop/README.md).


[TOC levels=1-3]: # ""

- [Aligning various genera from Protists](#aligning-various-genera-from-protists)
- [*Plasmodium*](#plasmodium)
  - [plasmodium: wgs](#plasmodium-wgs)
  - [plasmodium: assembly](#plasmodium-assembly)
  - [plasmodium: run](#plasmodium-run)


## Strain info

| Group             | Genus           | Genus ID | Comments | Species | Strains |
|:------------------|:----------------|---------:|:---------|--------:|--------:|
| Apicomplexans     |                 |          | 顶复虫    |         |         |
|                   | Plasmodium      |     5820 | 疟原虫属   |         |         |
|                   | Toxoplasma      |     5810 | 弓形虫属   |         |         |
|                   | Cryptosporidium |     5806 | 隐孢子虫属 |         |         |
|                   | Eimeria         |     5800 | 艾美球虫   |         |         |
|                   | Theileria       |     5873 | 泰勒虫属   |         |         |
|                   | Babesia         |     5864 | 巴倍虫属   |         |         |
| Oomycete          |                 |          | 卵菌      |         |         |
|                   | Phytophthora    |     4783 | 疫霉属    |         |         |
|                   | Pythium         |     4797 | 腐霉属    |         |         |
| Kinetoplastida    |                 |          | 动基体目   |         |         |
|                   | Leishmania      |     5658 | 利什曼虫属 |         |         |
|                   | Trypanosoma     |     5690 | 锥虫属    |         |         |
| Amoebozoa         |                 |          | 变形虫    |         |         |
|                   | Acanthamoeba    |     5754 | 棘阿米巴属 |         |         |
|                   | Entamoeba       |     5758 | 内阿米巴属 |         |         |
|                   | Dictyostelium   |     5782 | 网柄菌属   |         |         |
| Eustigmatophyceae |                 |          | 大眼藻纲   |         |         |
|                   | Nannochloropsis |     5748 | 微拟球藻   |         |         |
| Opalinata         |                 |          | 蛙片总纲   |         |         |
|                   | Blastocystis    |    12967 | 芽囊原虫属 |         |         |
| Metamonada        |                 |          | 后滴门    |         |         |
|                   | Giardia         |     5740 | 贾第虫属   |         |         |
| Euglenozoa        |                 |          | 眼虫门    |         |         |
|                   | Crithidia       |     5655 | 短膜虫属   |         |         |


## NCBI Assembly

```bash
export RANK_NAME=Protists

mkdir -p ~/data/alignment/${RANK_NAME}        # Working directory
cd ~/data/alignment/${RANK_NAME}

mysql -ualignDB -palignDB ar_refseq -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (
            5820, 5810, 5806, 5800, 5873,
            5864,
            4783, 4797,
            5658, 5690,
            5754, 5758, 5782, 
            5748, 
            12967,
            5740,
            5655
        )
    " \
    > raw.tsv

mysql -ualignDB -palignDB ar_genbank -e "
    SELECT 
        organism_name, species, genus, ftp_path, assembly_level
    FROM ar 
    WHERE 1=1
        AND genus_id in (
            5820, 5810, 5806, 5800, 5873,
            5864,
            4783, 4797,
            5658, 5690,
            5754, 5758, 5782, 
            5748, 
            12967,
            5740,
            5655
        )
    " \
    >> raw.tsv

cat raw.tsv |
    grep -v '^#' |
    tsv-filter --not-regex "1:\[.+\]" |
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
# vim ${RANK_NAME}.assembly.tsv
# cp ${RANK_NAME}.assembly.tsv ~/Scripts/withncbi/pop

# Cleaning
rm raw*.*sv

unset RANK_NAME

```

# *Plasmodium*


| name                      | taxon |
|:--------------------------|:------|
| Plasmodium                | 5820  |
| Plasmodium falciparum     | 5833  |
| Plasmodium falciparum 3D7 | 36329 |

Check NCBI pages

* http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=plasmodium&order=organism
* http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5820
* http://www.ncbi.nlm.nih.gov/assembly?term=txid5820[Organism:exp]
* http://www.ncbi.nlm.nih.gov/genome/?term=txid5820[Organism:exp]
* http://www.ncbi.nlm.nih.gov/genome/genomes/33

## plasmodium: wgs

Same as [here](README.md#wgstsv)

```bash
export RANK_LEVEL=genus
export RANK_ID=5820
export RANK_NAME=plasmodium

mkdir -p ~/data/alignment/${RANK_NAME}            # Working directory
cd ~/data/alignment/${RANK_NAME}

... # paste codes from README.md

# Cleaning
rm raw*.*sv

```

Edit them to fix names and comment out bad strains.

```bash
cd ~/data/alignment/plasmodium

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/withncbi/pop/plasmodium.wgs.tsv \
    --fix \
    -o WGS

bash WGS/plasmodium.wgs.rsync.sh

find WGS -name "*.gz" | parallel -j 4 gzip -t

```

## plasmodium: assembly

Same as [here](README.md#assemblytsv)

```bash
cd ~/data/alignment/${RANK_NAME}

... # paste codes from README.md

unset RANK_LEVEL
unset RANK_ID
unset RANK_NAME

```

```bash
cd ~/data/alignment/plasmodium

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/plasmodium.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/plasmodium.assembly.rsync.sh

bash ASSEMBLY/plasmodium.assembly.collect.sh

```

## plasmodium: run

```bash
$(brew --prefix)/Cellar/$(brew list --versions repeatmasker | sed 's/ /\//')/libexec/util/queryRepeatDatabase.pl \
    -species Alveolata -stat
```

* Rsync to hpcc

```bash
rsync -avP \
    ~/data/alignment/plasmodium/ \
    wangq@202.119.37.251:data/alignment/plasmodium

# rsync -avP wangq@202.119.37.251:data/alignment/plasmodium/ ~/data/alignment/plasmodium

```

`--perseq` for RefSeq Chromosome-level assemblies.

```bash
cd ~/data/alignment/plasmodium

# prep
egaz template \
    ASSEMBLY WGS \
    --prep -o GENOMES \
    --perseq Pfal_3D7 \
    --perseq Pber_ANKA \
    --perseq Pcha_chabaudi \
    --perseq Pcyn_strain_B \
    --perseq Pkno_strain_H \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Alveolata --parallel 24"

bsub -q mpi -n 24 -J "plasmodium-0_prep" "bash GENOMES/0_prep.sh"

ls -t output.* | head -n 1 | xargs tail -f | grep "==>"

# gff
for n in Pfal_3D7 Pber_ANKA Pcha_adami Pcha_chabaudi Pcyn_strain_B Pkno_strain_H; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"
    
    gzip -d -c ${FILE_GFF} > GENOMES/${n}/chr.gff
done

# multi
egaz template \
    GENOMES/Pfal_3D7 \
    $(find GENOMES -maxdepth 1 -type d -path "*/????*" | grep -v "Pfal_3D7") \
    --multi -o multi/ \
    --rawphylo --parallel 24 -v

bsub -q mpi -n 24 -J "plasmodium-1_pair" "bash multi/1_pair.sh"
bsub -w "ended(plasmodium-1_pair)" \
    -q mpi -n 24 -J "plasmodium-2_rawphylo" "bash multi/2_rawphylo.sh"
bsub  -w "ended(plasmodium-2_rawphylo)" \
    -q mpi -n 24 -J "plasmodium-3_multi" "bash multi/3_multi.sh"

# multi_Pfal
egaz template \
    GENOMES/Pfal_3D7 \
    $(find GENOMES -maxdepth 1 -type d -path "*/????*" | grep "Pfal_" | grep -v "Pfal_3D7") \
    GENOMES/Prei_SY57 \
    --multi -o multi/ \
    --multiname multi_Pfal --tree multi/Results/multi.nwk --outgroup Prei_SY57 \
    --parallel 24 -v

bsub -q mpi -n 24 -J "plasmodium-3_multi" "bash multi/3_multi.sh"

# self
egaz template \
    GENOMES/Pfal_3D7 GENOMES/Pber_ANKA GENOMES/Pcha_chabaudi \
    GENOMES/Pcyn_strain_B GENOMES/Pkno_strain_H \
    --self -o self/ \
    --circos --parallel 24 -v

bsub -q mpi -n 24 -J "plasmodium-1_self" "bash self/1_self.sh"
bsub -w "ended(plasmodium-1_self)" \
    -q mpi -n 24 -J "plasmodium-3_proc" "bash self/3_proc.sh"
bsub  -w "ended(plasmodium-3_proc)" \
    -q mpi -n 24 -J "plasmodium-4_circos" "bash self/4_circos.sh"

```

