# vsself alignDB for genomes from Ensembl

## `chr_length.csv`

```bash
cd ~/data/alignment/Ensembl

cat << 'EOF' > temp.tt
#!/bin/bash

cat << DELIMITER > chrUn.csv
common_name,taxon_id,chr,length,assembly
[% FOREACH item IN data -%]
[% item.name %],[% item.taxon %],chrUn,999999999,
[% END -%]
DELIMITER

if [ -f real_chr.csv ]; then
    rm real_chr.csv;
fi;

[% FOREACH item IN data -%]
# [% item.name %]
echo "==> [% item.name %]"
faops size [% item.name %]/*.fa > [% item.name %]/chr.sizes;
perl -aln -F"\t" -e 'print qq{[% item.name %],[% item.taxon %],$F[0],$F[1],}' [% item.name %]/chr.sizes >> real_chr.csv;

[% END -%]

cat chrUn.csv real_chr.csv > chr_length.csv
echo "==> chr_length.csv generated <=="

rm chrUn.csv
rm real_chr.csv
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
    --dir ~/data/alignment/Ensembl/Atha \
    --length 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl -d Athavsself -f ~/data/dumps/mysql/Athavsself.sql.gz
```

## OsatJap

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d OsatJapvsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d OsatJapvsself \
    -t OsatJap \
    --dir ~/data/alignment/Ensembl/OsatJap \
    --length 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl -d OsatJapvsself -f ~/data/dumps/mysql/OsatJapvsself.sql.gz
```

## Dmel

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d Dmelvsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d Dmelvsself \
    -t Dmel \
    --dir ~/data/alignment/Ensembl/Dmel \
    --length 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl -d Dmelvsself -f ~/data/dumps/mysql/Dmelvsself.sql.gz
```

## Mouse

```bash
perl ~/Scripts/alignDB/init/init_alignDB.pl \
    -d Mousevsself \
    --chr ~/data/alignment/Ensembl/chr_length.csv

perl ~/Scripts/alignDB/init/gen_alignDB_genome.pl \
    -d Mousevsself \
    -t Mouse \
    --dir ~/data/alignment/Ensembl/Mouse \
    --length 500000 \
    --parallel 8

perl ~/Scripts/alignDB/util/dup_db.pl -d Mousevsself -f ~/data/dumps/mysql/Mousevsself.sql.gz
```
