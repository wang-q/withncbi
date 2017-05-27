# `vsself` alignDB for genomes from Ensembl

[TOC levels=1-3]: # " "
- [`vsself` alignDB for genomes from Ensembl](#vsself-aligndb-for-genomes-from-ensembl)
    - [`chr_length.csv`](#chr-lengthcsv)
    - [Atha](#atha)
    - [OsatJap](#osatjap)
    - [Dmel](#dmel)
    - [Mouse](#mouse)
    - [S288c](#s288c)


## `chr_length.csv`

```bash
cd ~/data/alignment/Ensembl

cat << 'EOF' > temp.tt
#!/bin/bash

echo "common_name,taxon_id,chr,length,assembly" > chr_length.csv

[% FOREACH item IN data -%]
# [% item.name %]
echo "==> [% item.name %]"
faops size [% item.name %]/*.fa > [% item.name %]/chr.sizes;
perl -aln -F"\t" -e 'print qq{[% item.name %],[% item.taxon %],$F[0],$F[1],}' [% item.name %]/chr.sizes >> chr_length.csv;

[% END -%]

echo "==> chr_length.csv generated <=="

EOF

perl -MTemplate -e '
    my @data;
    for (grep {defined} @ARGV) {
        my ($n, $t) = split /=/;
        push @data, {name => $n, taxon => $t,};
    }
    my $tt = Template->new;
    $tt->process(
        q{temp.tt},
        {   data        => \@data,
            working_dir => $working_dir,
        },
    )
    ' \
    Cele=6239 Human=9606 OsatJap=39947 Atha=3702 Dmel=7227 Mouse=10090 S288c=559292 \
    > real_chr.sh

bash real_chr.sh
```

## Atha

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d Athavsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d Athavsself \
    -t Atha \
    --da ~/data/alignment/Ensembl/Atha \
    --truncate 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl --yes -d Athavsself -f ~/data/dumps/mysql/Athavsself.sql.gz
```

## OsatJap

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d OsatJapvsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d OsatJapvsself \
    -t OsatJap \
    --da ~/data/alignment/Ensembl/OsatJap \
    --truncate 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl --yes -d OsatJapvsself -f ~/data/dumps/mysql/OsatJapvsself.sql.gz
```

## Dmel

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d Dmelvsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d Dmelvsself \
    -t Dmel \
    --da ~/data/alignment/Ensembl/Dmel \
    --truncate 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl --yes -d Dmelvsself -f ~/data/dumps/mysql/Dmelvsself.sql.gz
```

## Mouse

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d Mousevsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d Mousevsself \
    -t Mouse \
    --da ~/data/alignment/Ensembl/Mouse \
    --truncate 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl --yes -d Mousevsself -f ~/data/dumps/mysql/Mousevsself.sql.gz
```

## S288c

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d S288cvsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d S288cvsself \
    -t S288c \
    --da ~/data/alignment/Ensembl/S288c \
    --truncate 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl --yes -d S288cvsself -f ~/data/dumps/mysql/S288cvsself.sql.gz
```
