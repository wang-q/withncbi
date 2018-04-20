# Aligning various genera from Protists

Less detailed than *Trichoderma* in
[README.md](https://github.com/wang-q/withncbi/blob/master/pop/README.md).


[TOC levels=1-3]: # " "
- [Aligning various genera from Protists](#aligning-various-genera-from-protists)
- [*Plasmodium*](#plasmodium)
    - [plasmodium: wgs](#plasmodium-wgs)
    - [plasmodium: assembly](#plasmodium-assembly)
    - [plasmodium: run](#plasmodium-run)


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

`--perseq` for Chromosome-level assemblies.

```bash
cd ~/data/alignment/plasmodium

# prep
egaz template \
    ASSEMBLY WGS \
    --prep -o GENOMES \
    --perseq Pfal_3D7 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Alveolata --parallel 24"

bsub -q mpi -n 24 -J "plasmodium-0_prep" "bash GENOMES/0_prep.sh"

ls -t output.* | head -n 1 | xargs tail -f | grep "==>"

# gff
for n in Pfal_3D7; do
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
    GENOMES/Calb_SC5314 \
    $(find GENOMES -maxdepth 1 -type d -path "*/????*" | grep "Pfal_" | grep -v "Pfal_3D7") \
    GENOMES/Cdub_CD36 \
    --multi -o multi/ \
    --multiname multi_Pfal --tree multi/Results/multi.nwk --outgroup Cdub_CD36 \
    --parallel 24 -v

bsub -q mpi -n 24 -J "plasmodium-3_multi" "bash multi/3_multi.sh"

# self
egaz template \
    GENOMES/Pfal_3D7 GENOMES/Cdub_CD36 GENOMES/Cort_Co_90_125  \
    --self -o self/ \
    --circos --parallel 24 -v

bsub -q mpi -n 24 -J "plasmodium-1_self" "bash self/1_self.sh"
bsub -w "ended(plasmodium-1_self)" \
    -q mpi -n 24 -J "plasmodium-3_proc" "bash self/3_proc.sh"
bsub  -w "ended(plasmodium-3_proc)" \
    -q mpi -n 24 -J "plasmodium-4_circos" "bash self/4_circos.sh"

```

