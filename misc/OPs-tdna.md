# T-DNA insertion sites

## Download

`Wed Jan 20 15:03:28 2016`

```bash
mkdir -p ~/data/salk
cd ~/data/salk

# Alternative url: http://signal.salk.edu/database/
perl ~/Scripts/download/list.pl -u http://natural.salk.edu/database/ -r transcriptome
perl ~/Scripts/download/download.pl -i database.yml -a
```

## Convert to positions

```bash
mkdir -p ~/data/salk/Atha
for name in CMT CSHL FLAG GABI IMAL MX RATM SAIL SALK SK WISC
do
    echo ${name};
    cat ~/data/salk/database/tdnaexpress/T-DNA.${name} \
        | perl -nla -e '
            @F >= 2 or next;
            next unless $F[1];
        
            my ( $chr, $pos ) = split /:/, $F[1];
            $chr =~ s/chr0?//i;
            $pos =~ s/^0+//;
            next unless $chr =~ /^\d+$/;
        
            print "$chr:$pos";
        ' \
        > ~/data/salk/Atha/T-DNA.${name}.pos.txt;
done

mkdir -p ~/data/salk/OstaJap
for name in affjp cirad csiro genoplante gsnu ostid PFG_FSTs rifgp rmd ship trim ucd
do
    echo ${name};
    cat ~/data/salk/database/RiceGE/T-DNA.${name} \
        | perl -nla -e '
            @F >= 2 or next;
            next unless $F[1];
        
            my ( $chr, $pos ) = split /:/, $F[1];
            $chr =~ s/chr0?//i;
            $pos =~ s/^0+//;
            next unless $chr =~ /^\d+$/;
        
            print "$chr:$pos";
        ' \
        > ~/data/salk/OstaJap/T-DNA.${name}.pos.txt;
done

```

## Restore vsself and pop databases

[`vsself`](../pop/OPs-vsself.md) and [`pop`](../pop/OPs-multi.md)

```bash
perl ~/Scripts/alignDB/util/dup_db.pl -g Athavsself_tdna -f ~/data/dumps/mysql/Athavsself.sql.gz
perl ~/Scripts/alignDB/util/dup_db.pl -g Ath_n19_pop_tdna -f ~/data/dumps/mysql/Ath_n19_pop.sql.gz
```

## Use style 'center'

```bash
cd ~/data/salk/Atha

perl ~/Scripts/alignDB/init/insert_position.pl \
    -d Athavsself_tdna  --style center \
    --batch 50 --parallel 8 \
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

perl ~/Scripts/alignDB/init/update_sw_cv.pl \
    -d Athavsself_tdna \
    --batch 10 --parallel 8
perl ~/Scripts/alignDB/init/update_annotation.pl \
    -d Athavsself_tdna \
    -a ~/data/alignment/Ensembl/Atha/anno.yml \
    --batch 10 --parallel 8

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
    --by tt --index --chart \
    -d Athavsself_tdna -o ~/data/salk/Athavsself_tdna.ofg.xlsx
```
