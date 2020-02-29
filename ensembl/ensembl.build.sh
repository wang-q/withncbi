#!/usr/bin/env bash

# Arabidopsis thaliana
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db arabidopsis_thaliana_core_45_98_11 \
     --ensembl /home/wangq/data/ensembl98/mysql/arabidopsis_thaliana_core_45_98_11

# Aspergillus fumigatus
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db aspergillus_fumigatus_core_45_98_1 \
     --ensembl /home/wangq/data/ensembl98/mysql/aspergillus_fumigatus_core_45_98_1

# Caenorhabditis elegans
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db caenorhabditis_elegans_core_98_269 \
     --ensembl /home/wangq/data/ensembl98/mysql/caenorhabditis_elegans_core_98_269

# Dictyostelium discoideum
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db dictyostelium_discoideum_core_45_98_1 \
     --ensembl /home/wangq/data/ensembl98/mysql/dictyostelium_discoideum_core_45_98_1

# Drosophila melanogaster
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db drosophila_melanogaster_core_98_7 \
     --ensembl /home/wangq/data/ensembl98/mysql/drosophila_melanogaster_core_98_7

# Homo sapiens
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db homo_sapiens_core_98_38 \
     --ensembl /home/wangq/data/ensembl98/mysql/homo_sapiens_core_98_38

# Mus musculus
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db mus_musculus_core_98_38 \
     --ensembl /home/wangq/data/ensembl98/mysql/mus_musculus_core_98_38

# Oryza sativa
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db oryza_sativa_core_45_98_7 \
     --ensembl /home/wangq/data/ensembl98/mysql/oryza_sativa_core_45_98_7

# Plasmodium falciparum
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db plasmodium_falciparum_core_45_98_1 \
     --ensembl /home/wangq/data/ensembl98/mysql/plasmodium_falciparum_core_45_98_1

# Saccharomyces cerevisiae
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db saccharomyces_cerevisiae_core_98_4 \
     --ensembl /home/wangq/data/ensembl98/mysql/saccharomyces_cerevisiae_core_98_4

# Schizosaccharomyces pombe
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db schizosaccharomyces_pombe_core_45_98_2 \
     --ensembl /home/wangq/data/ensembl98/mysql/schizosaccharomyces_pombe_core_45_98_2

