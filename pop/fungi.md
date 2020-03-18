# Aligning various genera from Fungi

Less detailed than *Trichoderma* in
[README.md](https://github.com/wang-q/withncbi/blob/master/pop/README.md).


[TOC levels=1-3]: # ""

- [Aligning various genera from Fungi](#aligning-various-genera-from-fungi)
- [*Candida*](#candida)
  - [candida: wgs](#candida-wgs)
  - [candida: assembly](#candida-assembly)
  - [candida: phylo](#candida-phylo)
  - [candida: run](#candida-run)
- [*Penicillium*](#penicillium)
  - [penicillium: wgs](#penicillium-wgs)
  - [penicillium: assembly](#penicillium-assembly)
  - [penicillium: phylo](#penicillium-phylo)


# *Candida*

Check NCBI pages

* https://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Candida*&order=organism
* https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1535326
* https://www.ncbi.nlm.nih.gov/assembly?term=txid1535326[Organism:exp]
* https://www.ncbi.nlm.nih.gov/genome/?term=txid1535326[Organism:exp]

Candida dubliniensis CD36 (11x sanger)
Candida albicans WO-1 as target (10x sanger)
Candida tropicalis MYA-3404 (8.85x sanger)

## candida: wgs

Same as [here](README.md#wgstsv)

```bash
export RANK_LEVEL=genus
export RANK_ID=1535326
export RANK_NAME=candida

mkdir -p ~/data/alignment/${RANK_NAME}            # Working directory
cd ~/data/alignment/${RANK_NAME}

... # paste codes from README.md

# Cleaning
rm raw*.*sv

```

Edit them to fix names and comment out bad strains.

```bash
cd ~/data/alignment/candida

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/withncbi/pop/candida.wgs.tsv \
    --fix \
    -o WGS

bash WGS/candida.wgs.rsync.sh
bash WGS/candida.wgs.aria2.sh

find WGS -name "*.gz" | parallel -j 4 gzip -t

```

## candida: assembly

Same as [here](README.md#assemblytsv)

```bash
cd ~/data/alignment/${RANK_NAME}

... # paste codes from README.md

unset RANK_LEVEL
unset RANK_ID
unset RANK_NAME

```

```bash
cd ~/data/alignment/candida

perl ~/Scripts/withncbi/taxon/assembly_prep.pl \
    -f ~/Scripts/withncbi/pop/candida.assembly.tsv \
    -o ASSEMBLY

bash ASSEMBLY/candida.assembly.rsync.sh

bash ASSEMBLY/candida.assembly.collect.sh

```

## candida: phylo

Same as [here](README.md#mash)

```bash
mkdir -p ~/data/alignment/candida/mash
cd ~/data/alignment/candida/mash

... # paste codes from README.md


nw_display -w 600 -b 'visibility:hidden' -s tree.nwk |
    rsvg-convert -o ~/Scripts/withncbi/image/Candida.png

```

![Candida.png](../image/Candida.png)

## candida: run

* Rsync to hpcc

```bash
rsync -avP \
    ~/data/alignment/candida/ \
    wangq@202.119.37.251:data/alignment/candida

# rsync -avP wangq@202.119.37.251:data/alignment/candida/ ~/data/alignment/candida

```

`--perseq` for Chromosome-level assemblies.

```bash
cd ~/data/alignment/candida

# prep
egaz template \
    ASSEMBLY WGS \
    --prep -o GENOMES \
    --perseq Calb_WO_1 \
    --perseq Cdub_CD36 \
    --perseq Ctro_MYA_3404 \
    --min 5000 --about 5000000 \
    -v --repeatmasker "--species Fungi --parallel 24"

bsub -q mpi -n 24 -J "candida-0_prep" "bash GENOMES/0_prep.sh"

ls -t output.* | head -n 1 | xargs tail -f | grep "==>"

# gff
for n in Calb_WO_1 Cdub_CD36 Ctro_MYA_3404; do
    FILE_GFF=$(find ASSEMBLY -type f -name "*_genomic.gff.gz" | grep "${n}")
    echo >&2 "==> Processing ${n}/${FILE_GFF}"
    
    gzip -d -c ${FILE_GFF} > GENOMES/${n}/chr.gff
done

# sanger
egaz template \
    GENOMES/Calb_WO_1 \
    GENOMES/Cdub_CD36 \
    GENOMES/Ctro_MYA_3404 \
    --multi -o multi/ \
    --multiname sanger --order \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_pair" "bash multi/1_pair.sh"
bsub  -w "ended(candida-1_pair)" \
    -q mpi -n 24 -J "candida-3_multi" "bash multi/3_multi.sh"

# multi
egaz template \
    GENOMES/Calb_WO_1 \
    $(find GENOMES -maxdepth 1 -mindepth 1 -type d | grep -v "Calb_WO_1") \
    --multi -o multi/ \
    --tree mash/tree.nwk \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_pair" "bash multi/1_pair.sh"
bsub  -w "ended(candida-1_pair)" \
    -q mpi -n 24 -J "candida-3_multi" "bash multi/3_multi.sh"

# multi_Calb
egaz template \
    GENOMES/Calb_WO_1 \
    $(find GENOMES -maxdepth 1 -mindepth 1 -type d | grep "Calb_" | grep -v "Calb_SC5314") \
    GENOMES/Cdub_CD36 \
    --multi -o multi/ \
    --multiname Calb --tree mash/tree.nwk --outgroup Cdub_CD36 \
    --parallel 24 -v

bsub -q mpi -n 24 -J "candida-3_multi" "bash multi/3_multi.sh"

# self
egaz template \
    GENOMES/Calb_WO_1 \
    GENOMES/Cdub_CD36 \
    GENOMES/Ctro_MYA_3404 \
    --self -o self/ \
    --circos --parallel 24 -v

bsub -q mpi -n 24 -J "candida-1_self" "bash self/1_self.sh"
bsub -w "ended(candida-1_self)" \
    -q mpi -n 24 -J "candida-3_proc" "bash self/3_proc.sh"
bsub  -w "ended(candida-3_proc)" \
    -q mpi -n 24 -J "candida-4_circos" "bash self/4_circos.sh"

```

# *Penicillium*

* http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=penicillium&order=organism
* http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5073
* http://www.ncbi.nlm.nih.gov/assembly?term=txid5073[Organism:exp]
* http://www.ncbi.nlm.nih.gov/genome/?term=txid5073[Organism:exp]

## penicillium: wgs

Same as [here](README.md#poptrichodermawgstsv)

```bash
export RANK_LEVEL=genus
export RANK_ID=5073
export RANK_NAME=penicillium

mkdir -p ~/data/alignment/${RANK_NAME}            # Working directory
cd ~/data/alignment/${RANK_NAME}

... # paste codes from README.md

# Cleaning
rm raw*.*sv

```

Edit them to fix names and comment out bad strains.

```bash
cd ~/data/alignment/penicillium

perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
    -f ~/Scripts/withncbi/pop/penicillium.wgs.tsv \
    --fix \
    -o WGS \
    -a

aria2c -UWget -x 6 -s 3 -c -i WGS/penicillium.wgs.url.txt

find WGS -name "*.gz" | xargs gzip -t

```

## penicillium: assembly

Same as [here](README.md#assembly_preppl)

```bash
cd ~/data/alignment/${RANK_NAME}

... # paste codes from README.md

unset RANK_LEVEL
unset RANK_ID
unset RANK_NAME

```

```bash
cd ~/data/alignment/penicillium

bash ASSEMBLY/penicillium.assembly.rsync.sh

bash ASSEMBLY/penicillium.assembly.collect.sh

```

There're no good target. Pchr_P2niaD18 is the only one on chromosome level, but is not de novo
assembled and hasn't annotations.
