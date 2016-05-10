# T-DNA insertion sites

## Download

Wed Jan 20 15:03:28 2016

```bash
mkdir -p ~/data/salk
cd ~/data/salk

# Alternative url: http://signal.salk.edu/database/
perl ~/Scripts/download/list.pl -u http://natural.salk.edu/database/ -r transcriptome
perl ~/Scripts/download/download.pl -i database.yml -a
```

## Convert to beds

```bash
mkdir -p ~/data/salk/Atha
mkdir -p ~/data/salk/OstaJap

for name in CMT CSHL FLAG GABI IMAL MX RATM SAIL SALK SK WISC
do
    echo $name;
    perl ~/Scripts/alignDB/ofg/tdna2position.pl -i ~/data/salk/database/tdnaexpress/T-DNA.$name -o ~/data/salk/Atha/T-DNA.$name.pos.txt;
done

for name in affjp cirad csiro genoplante gsnu ostid PFG_FSTs rifgp rmd ship trim ucd
do
    echo $name;
    perl ~/Scripts/alignDB/ofg/tdna2position.pl -i ~/data/salk/database/RiceGE/T-DNA.$name -o ~/data/salk/OstaJap/T-DNA.$name.pos.txt;
done

```

## Restore databases

```bash
perl ~/Scripts/alignDB/util/dup_db.pl -g Athavsself_tdna -f ~/data/dumps/mysql/Athavsself.sql.gz

```

## Use style 'center'

```bash
perl ~/Scripts/alignDB/ofg/insert_position.pl \
    -d Athavsself_tdna  --style center --batch 50 --parallel 8 \
    --tag tdna --type CMT  -f ~/data/salk/Atha/T-DNA.CMT.pos.txt  \
    --tag tdna --type CSHL -f ~/data/salk/Atha/T-DNA.CSHL.pos.txt \
    --tag tdna --type FLAG -f ~/data/salk/Atha/T-DNA.FLAG.pos.txt \
    --tag tdna --type GABI -f ~/data/salk/Atha/T-DNA.GABI.pos.txt \
    --tag tdna --type IMAL -f ~/data/salk/Atha/T-DNA.IMAL.pos.txt \
    --tag tdna --type MX   -f ~/data/salk/Atha/T-DNA.MX.pos.txt   \
    --tag tdna --type RATM -f ~/data/salk/Atha/T-DNA.RATM.pos.txt \
    --tag tdna --type SAIL -f ~/data/salk/Atha/T-DNA.SAIL.pos.txt \
    --tag tdna --type SALK -f ~/data/salk/Atha/T-DNA.SALK.pos.txt \
    --tag tdna --type SK   -f ~/data/salk/Atha/T-DNA.SK.pos.txt   \
    --tag tdna --type WISC -f ~/data/salk/Atha/T-DNA.WISC.pos.txt

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d Athavsself_tdna --batch 10 --parallel 8
perl ~/Scripts/alignDB/init/update_feature.pl -d Athavsself_tdna \
    -e arabidopsis_thaliana_core_29_82_10 --batch 10 --parallel 8

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl --by tt -d Athavsself_tdna -o Athavsself_tdna.ofg.xlsx

```
