# Aligning various genera from Fungi

Less detailed than *Trichoderma* in
[README.md](https://github.com/wang-q/withncbi/blob/master/pop/README.md).


[TOC levels=1-3]: # " "
- [Aligning various genera from Fungi](#aligning-various-genera-from-fungi)
- [*Candida*](#candida)
    - [candida: wgs](#candida-wgs)
    - [candida: assembly](#candida-assembly)
- [*Penicillium*](#penicillium)
    - [penicillium: wgs](#penicillium-wgs)
    - [penicillium: assembly](#penicillium-assembly)


# *Candida*

Check NCBI pages

* http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Candida*&order=organism
* http://www.ncbi.nlm.nih.gov/assembly?term=txid1535326[Organism:exp]
* http://www.ncbi.nlm.nih.gov/genome/?term=txid1535326[Organism:exp]

## candida: wgs

Same as [here](README.md#poptrichodermawgstsv)

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

find WGS -name "*.gz" | parallel -j 4 gzip -t

```


## candida: assembly

Same as [here](README.md#assembly_preppl)

```bash
cd ~/data/alignment/${RANK_NAME}

... # paste codes from README.md

unset RANK_LEVEL
unset RANK_ID
unset RANK_NAME

```

```bash
cd ~/data/alignment/candida

bash ASSEMBLY/candida.assembly.rsync.sh

# linode
# rsync -avP wangq@173.230.144.105:data/alignment/candida/ ~/data/alignment/candida

bash ASSEMBLY/candida.assembly.collect.sh

```

    --plan 'name=four_way;t=Cdub_CD36;qs=Corh_Co_90_125,Calb_WO_1,Ctro_MYA_3404'
    --plan 'name=four_way_2;t=Corh_Co_90_125;qs=Cdub_CD36,Calb_WO_1,Ctro_MYA_3404'

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

# rsync -avP wangq@173.230.144.105:data/alignment/penicillium/ ~/data/alignment/penicillium

bash ASSEMBLY/penicillium.assembly.collect.sh

```

There're no good target. Pchr_P2niaD18 is the only one on chromosome level, but is not de novo
assembled and hasn't annotations.
