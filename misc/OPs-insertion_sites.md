# Transposons insertion sites

## Fly ([the Gene Disruption Project](http://flypush.imgen.bcm.tmc.edu/pscreen/))

### Get data

Original paper:

MI, RMCE collection (by injection):
Nagarkar-Jaiswal S, Lee PT, Campbell ME, Chen K, Anguiano-Zarate S, Cantu Gutierrez M, Busby T, Lin
WW, He Y, Schulze KL, Booth BW, Evans-Holm M, Venken KJ, Levis RW, Spradling AC, Hoskins RA, Bellen
HJ (2015) A library of MiMICs allows tagging of genes and reversible spatial and temporal knockdown
of proteins in Drosophila. eLife 4:e05338.

MI, RMCE technique:
Venken KJT, Schulze KL, Haelterman NA, Pan H, He Y, Evans-Holm M, Carlson JW, Levis RW, Spradling
AC, Hoskins RA, Bellen HJ (2011) MiMIC: a highly versatile transposon insertion resource for
engineering Drosophila melanogaster genes. Nature Methods 8:737-743.


MB, EY and piggyBac:
Bellen HJ, Levis RW, He Y, Carlson JW, Evans-Holm M, Bae E, Kim J, Metaxakis A, Savakis C, Schulze
KL, Hoskins RA, Spradling AC (2011) The Drosophila Gene Disruption Project: progress using
transposons with distinctive site-specificities. Genetics 188:731-743.

```bash
mkdir -p ~/data/ofg/gdp
cd ~/data/ofg/gdp

wget -N http://flypush.imgen.bcm.tmc.edu/pscreen/files/MI-list-2015-07-27.xlsx

perl ~/Scripts/fig_table/xlsx2csv.pl -f MI-list-2015-07-27.xlsx --sheet 'Sheet1' \
    | perl -F/,/ -an -e '
        /^Insertion/i and next;
        my $loc = $F[3];
        $loc =~ s/"//g;
        $loc =~ s/\s*\[.*\]//;
        if ($loc =~ /^\w+:\d+$/) {
            print qq{$loc\n};
        }
    ' \
    | sort \
    > MI.pos.txt

wget -N http://flypush.imgen.bcm.tmc.edu/pscreen/files/MB_sites_Bellen_2011.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f MB_sites_Bellen_2011.xls \
    | perl -F/,/ -an -e '
        /^Line/i and next;
        my $chr = $F[2];
        my $loc = $F[3];
        print qq{$chr:$loc\n};
    ' \
    | sort \
    > MB.pos.txt

wget -N http://flypush.imgen.bcm.tmc.edu/pscreen/files/EY_sites_Bellen_2011.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f EY_sites_Bellen_2011.xls \
    | perl -F/,/ -an -e '
        /^Line/i and next;
        my $chr = $F[2];
        my $loc = $F[3];
        print qq{$chr:$loc\n};
    ' \
    | sort \
    > EY.pos.txt

wget -N http://flypush.imgen.bcm.tmc.edu/pscreen/files/Exelixis_PBac_sites_Bellen_2011.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f Exelixis_PBac_sites_Bellen_2011.xls \
    | perl -F/,/ -an -e '
        /^Line/i and next;
        my $chr = $F[2];
        my $loc = $F[3];
        print qq{$chr:$loc\n};
    ' \
    | sort \
    > PiggyBac.pos.txt

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
