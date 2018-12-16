#!/usr/bin/env bash

# Arabidopsis thaliana
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db arabidopsis_thaliana_core_41_94_11 \
     --ensembl /home/wangq/data/ensembl94/mysql/arabidopsis_thaliana_core_41_94_11

# Aspergillus fumigatus
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db aspergillus_fumigatus_core_41_94_1 \
     --ensembl /home/wangq/data/ensembl94/mysql/aspergillus_fumigatus_core_41_94_1

# Caenorhabditis elegans
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db caenorhabditis_elegans_core_94_260 \
     --ensembl /home/wangq/data/ensembl94/mysql/caenorhabditis_elegans_core_94_260

# Dictyostelium discoideum
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db dictyostelium_discoideum_core_41_94_1 \
     --ensembl /home/wangq/data/ensembl94/mysql/dictyostelium_discoideum_core_41_94_1

# Drosophila melanogaster
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db drosophila_melanogaster_core_94_6 \
     --ensembl /home/wangq/data/ensembl94/mysql/drosophila_melanogaster_core_94_6

# Homo sapiens
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db homo_sapiens_core_94_37 \
     --ensembl /home/wangq/data/ensembl94/mysql/homo_sapiens_core_94_37

# Mus musculus
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db mus_musculus_core_94_38 \
     --ensembl /home/wangq/data/ensembl94/mysql/mus_musculus_core_94_38

# Oryza sativa
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db oryza_sativa_core_41_94_7 \
     --ensembl /home/wangq/data/ensembl94/mysql/oryza_sativa_core_41_94_7

# Plasmodium falciparum
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db plasmodium_falciparum_core_41_94_1 \
     --ensembl /home/wangq/data/ensembl94/mysql/plasmodium_falciparum_core_41_94_1

# Saccharomyces cerevisiae
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db saccharomyces_cerevisiae_core_94_4 \
     --ensembl /home/wangq/data/ensembl94/mysql/saccharomyces_cerevisiae_core_94_4

# Schizosaccharomyces pombe
perl /home/wangq/Scripts/withncbi/ensembl/build_ensembl.pl \
     -s localhost --port 3306 -u alignDB --password alignDB \
     --initdb --db schizosaccharomyces_pombe_core_41_94_2 \
     --ensembl /home/wangq/data/ensembl94/mysql/schizosaccharomyces_pombe_core_41_94_2

