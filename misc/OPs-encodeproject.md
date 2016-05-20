# [Encode Project](https://www.encodeproject.org/)

## Data with narrowPeak

Data of this type are naturally fit for OFG.

```bash
mkdir -p ~/data/encodeproject
cd ~/data/encodeproject

curl -o report.narrowPeak.tsv 'https://www.encodeproject.org/report.tsv?type=Experiment&files.file_type=bed+narrowPeak&status=released'
curl -o metadata.narrowPeak.tsv 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=bed+narrowPeak&status=released/metadata.tsv'

```


### Process

1. Dmelvsself center

    ```bash
    cd ~/data/ofg/gdp

    perl ~/Scripts/alignDB/util/dup_db.pl -g Dmelvsself_gdp -f ~/data/dumps/mysql/Dmelvsself.sql.gz

    perl ~/Scripts/alignDB/ofg/insert_position.pl \
        -d Dmelvsself_gdp  --style center --batch 50 --parallel 8 \
        --tag transposon --type MI       -f ~/data/ofg/gdp/MI.pos.txt \
        --tag transposon --type MB       -f ~/data/ofg/gdp/MB.pos.txt \
        --tag transposon --type EY       -f ~/data/ofg/gdp/EY.pos.txt \
        --tag transposon --type PiggyBac -f ~/data/ofg/gdp/PiggyBac.pos.txt

    perl ~/Scripts/alignDB/init/update_sw_cv.pl \
        -d Dmelvsself_gdp \
        --batch 10 --parallel 8
    perl ~/Scripts/alignDB/init/update_annotation.pl \
        -d Dmelvsself_gdp \
        -a ~/data/alignment/Ensembl/Dmel/anno.yml \
        --batch 10 --parallel 8

    perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
        --by tt --index --chart \
        --replace ofg="insertion sites" \
        -d Dmelvsself_gdp -o ~/data/salk/Dmelvsself_gdp.ofg.xlsx
    ```

2. Dmel_n22_pop center

    ```bash
    cd ~/data/ofg/gdp

    perl ~/Scripts/alignDB/util/dup_db.pl -g Dmel_n22_pop_gdp -f ~/data/dumps/mysql/Dmel_n22_pop.sql.gz

    perl ~/Scripts/alignDB/ofg/insert_position.pl \
        -d Dmel_n22_pop_gdp  --style center --batch 50 --parallel 8 \
        --tag transposon --type MI       -f ~/data/ofg/gdp/MI.pos.txt \
        --tag transposon --type MB       -f ~/data/ofg/gdp/MB.pos.txt \
        --tag transposon --type EY       -f ~/data/ofg/gdp/EY.pos.txt \
        --tag transposon --type PiggyBac -f ~/data/ofg/gdp/PiggyBac.pos.txt

    perl ~/Scripts/alignDB/init/update_sw_cv.pl \
        -d Dmel_n22_pop_gdp \
        --batch 10 --parallel 8
    perl ~/Scripts/alignDB/init/update_annotation.pl \
        -d Dmel_n22_pop_gdp \
        -a ~/data/alignment/Ensembl/Dmel/anno.yml \
        --batch 10 --parallel 8

    perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl \
        --by tt --index --chart \
        --replace ofg="insertion sites" \
        -d Dmel_n22_pop_gdp -o ~/data/salk/Dmel_n22_pop_gdp.ofg.xlsx
    ```
