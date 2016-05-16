#!/usr/bin/env bash

# Arabidopsis thaliana
if [ -d /home/wangq/data/alignment/Ensembl/Atha ];
then
    echo "==> Arabidopsis thaliana"

    if [ -f /home/wangq/data/alignment/Ensembl/Atha/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Atha/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use arabidopsis_thaliana_core_29_82_10'; then
        cd /home/wangq/data/alignment/Ensembl/Atha

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e arabidopsis_thaliana_core_29_82_10 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e arabidopsis_thaliana_core_29_82_10 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database arabidopsis_thaliana_core_29_82_10 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Atha does not exist"
fi

# Aspergillus fumigatus
if [ -d /home/wangq/data/alignment/Ensembl/Afum ];
then
    echo "==> Aspergillus fumigatus"

    if [ -f /home/wangq/data/alignment/Ensembl/Afum/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Afum/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use aspergillus_fumigatus_core_29_82_2'; then
        cd /home/wangq/data/alignment/Ensembl/Afum

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e aspergillus_fumigatus_core_29_82_2 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e aspergillus_fumigatus_core_29_82_2 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database aspergillus_fumigatus_core_29_82_2 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Afum does not exist"
fi

# Brassica rapa
if [ -d /home/wangq/data/alignment/Ensembl/Brap ];
then
    echo "==> Brassica rapa"

    if [ -f /home/wangq/data/alignment/Ensembl/Brap/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Brap/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use brassica_rapa_core_29_82_1'; then
        cd /home/wangq/data/alignment/Ensembl/Brap

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e brassica_rapa_core_29_82_1 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e brassica_rapa_core_29_82_1 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database brassica_rapa_core_29_82_1 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Brap does not exist"
fi

# Caenorhabditis elegans
if [ -d /home/wangq/data/alignment/Ensembl/Cele ];
then
    echo "==> Caenorhabditis elegans"

    if [ -f /home/wangq/data/alignment/Ensembl/Cele/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Cele/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use caenorhabditis_elegans_core_82_245'; then
        cd /home/wangq/data/alignment/Ensembl/Cele

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e caenorhabditis_elegans_core_82_245 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e caenorhabditis_elegans_core_82_245 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database caenorhabditis_elegans_core_82_245 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Cele does not exist"
fi

# Dictyostelium discoideum
if [ -d /home/wangq/data/alignment/Ensembl/Ddis ];
then
    echo "==> Dictyostelium discoideum"

    if [ -f /home/wangq/data/alignment/Ensembl/Ddis/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Ddis/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use dictyostelium_discoideum_core_29_82_1'; then
        cd /home/wangq/data/alignment/Ensembl/Ddis

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e dictyostelium_discoideum_core_29_82_1 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e dictyostelium_discoideum_core_29_82_1 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database dictyostelium_discoideum_core_29_82_1 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Ddis does not exist"
fi

# Drosophila melanogaster
if [ -d /home/wangq/data/alignment/Ensembl/Dmel ];
then
    echo "==> Drosophila melanogaster"

    if [ -f /home/wangq/data/alignment/Ensembl/Dmel/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Dmel/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use drosophila_melanogaster_core_82_602'; then
        cd /home/wangq/data/alignment/Ensembl/Dmel

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e drosophila_melanogaster_core_82_602 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e drosophila_melanogaster_core_82_602 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database drosophila_melanogaster_core_82_602 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Dmel does not exist"
fi

# Drosophila simulans
if [ -d /home/wangq/data/alignment/Ensembl/Dsim ];
then
    echo "==> Drosophila simulans"

    if [ -f /home/wangq/data/alignment/Ensembl/Dsim/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Dsim/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use drosophila_simulans_core_29_82_1'; then
        cd /home/wangq/data/alignment/Ensembl/Dsim

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e drosophila_simulans_core_29_82_1 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e drosophila_simulans_core_29_82_1 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database drosophila_simulans_core_29_82_1 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Dsim does not exist"
fi

# Gorilla gorilla
if [ -d /home/wangq/data/alignment/Ensembl/Gorilla ];
then
    echo "==> Gorilla gorilla"

    if [ -f /home/wangq/data/alignment/Ensembl/Gorilla/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Gorilla/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use gorilla_gorilla_core_82_31'; then
        cd /home/wangq/data/alignment/Ensembl/Gorilla

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e gorilla_gorilla_core_82_31 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e gorilla_gorilla_core_82_31 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database gorilla_gorilla_core_82_31 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Gorilla does not exist"
fi

# Homo sapiens
if [ -d /home/wangq/data/alignment/Ensembl/Human ];
then
    echo "==> Homo sapiens"

    if [ -f /home/wangq/data/alignment/Ensembl/Human/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Human/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use homo_sapiens_core_82_37'; then
        cd /home/wangq/data/alignment/Ensembl/Human

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e homo_sapiens_core_82_37 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e homo_sapiens_core_82_37 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database homo_sapiens_core_82_37 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Human does not exist"
fi

# Macaca mulatta
if [ -d /home/wangq/data/alignment/Ensembl/Rhesus ];
then
    echo "==> Macaca mulatta"

    if [ -f /home/wangq/data/alignment/Ensembl/Rhesus/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Rhesus/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use macaca_mulatta_core_82_10'; then
        cd /home/wangq/data/alignment/Ensembl/Rhesus

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e macaca_mulatta_core_82_10 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e macaca_mulatta_core_82_10 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database macaca_mulatta_core_82_10 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Rhesus does not exist"
fi

# Mus musculus
if [ -d /home/wangq/data/alignment/Ensembl/Mouse ];
then
    echo "==> Mus musculus"

    if [ -f /home/wangq/data/alignment/Ensembl/Mouse/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Mouse/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use mus_musculus_core_82_38'; then
        cd /home/wangq/data/alignment/Ensembl/Mouse

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e mus_musculus_core_82_38 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e mus_musculus_core_82_38 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database mus_musculus_core_82_38 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Mouse does not exist"
fi

# Oryza sativa
if [ -d /home/wangq/data/alignment/Ensembl/OsatJap ];
then
    echo "==> Oryza sativa"

    if [ -f /home/wangq/data/alignment/Ensembl/OsatJap/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/OsatJap/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use oryza_sativa_core_29_82_7'; then
        cd /home/wangq/data/alignment/Ensembl/OsatJap

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e oryza_sativa_core_29_82_7 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e oryza_sativa_core_29_82_7 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database oryza_sativa_core_29_82_7 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatJap does not exist"
fi

# Pan troglodytes
if [ -d /home/wangq/data/alignment/Ensembl/Chimp ];
then
    echo "==> Pan troglodytes"

    if [ -f /home/wangq/data/alignment/Ensembl/Chimp/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Chimp/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use pan_troglodytes_core_82_214'; then
        cd /home/wangq/data/alignment/Ensembl/Chimp

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e pan_troglodytes_core_82_214 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e pan_troglodytes_core_82_214 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database pan_troglodytes_core_82_214 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Chimp does not exist"
fi

# Plasmodium falciparum
if [ -d /home/wangq/data/alignment/Ensembl/Pfal ];
then
    echo "==> Plasmodium falciparum"

    if [ -f /home/wangq/data/alignment/Ensembl/Pfal/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Pfal/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use plasmodium_falciparum_core_29_82_3'; then
        cd /home/wangq/data/alignment/Ensembl/Pfal

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e plasmodium_falciparum_core_29_82_3 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e plasmodium_falciparum_core_29_82_3 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database plasmodium_falciparum_core_29_82_3 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Pfal does not exist"
fi

# Pongo abelii
if [ -d /home/wangq/data/alignment/Ensembl/Orangutan ];
then
    echo "==> Pongo abelii"

    if [ -f /home/wangq/data/alignment/Ensembl/Orangutan/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Orangutan/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use pongo_abelii_core_82_1'; then
        cd /home/wangq/data/alignment/Ensembl/Orangutan

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e pongo_abelii_core_82_1 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e pongo_abelii_core_82_1 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database pongo_abelii_core_82_1 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Orangutan does not exist"
fi

# Saccharomyces cerevisiae
if [ -d /home/wangq/data/alignment/Ensembl/S288c ];
then
    echo "==> Saccharomyces cerevisiae"

    if [ -f /home/wangq/data/alignment/Ensembl/S288c/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/S288c/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use saccharomyces_cerevisiae_core_29_82_4'; then
        cd /home/wangq/data/alignment/Ensembl/S288c

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e saccharomyces_cerevisiae_core_29_82_4 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e saccharomyces_cerevisiae_core_29_82_4 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database saccharomyces_cerevisiae_core_29_82_4 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/S288c does not exist"
fi

# Schizosaccharomyces pombe
if [ -d /home/wangq/data/alignment/Ensembl/Spom ];
then
    echo "==> Schizosaccharomyces pombe"

    if [ -f /home/wangq/data/alignment/Ensembl/Spom/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Spom/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use schizosaccharomyces_pombe_core_29_82_2'; then
        cd /home/wangq/data/alignment/Ensembl/Spom

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e schizosaccharomyces_pombe_core_29_82_2 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e schizosaccharomyces_pombe_core_29_82_2 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database schizosaccharomyces_pombe_core_29_82_2 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Spom does not exist"
fi

# Solanum lycopersicum
if [ -d /home/wangq/data/alignment/Ensembl/Slyc ];
then
    echo "==> Solanum lycopersicum"

    if [ -f /home/wangq/data/alignment/Ensembl/Slyc/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Slyc/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use solanum_lycopersicum_core_29_82_250'; then
        cd /home/wangq/data/alignment/Ensembl/Slyc

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e solanum_lycopersicum_core_29_82_250 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e solanum_lycopersicum_core_29_82_250 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database solanum_lycopersicum_core_29_82_250 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Slyc does not exist"
fi

# Solanum tuberosum
if [ -d /home/wangq/data/alignment/Ensembl/Stub ];
then
    echo "==> Solanum tuberosum"

    if [ -f /home/wangq/data/alignment/Ensembl/Stub/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Stub/anno.yml exists"
    elif mysql -hlocalhost -ualignDB -palignDB -e 'use solanum_tuberosum_core_29_82_4'; then
        cd /home/wangq/data/alignment/Ensembl/Stub

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e solanum_tuberosum_core_29_82_4 \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e solanum_tuberosum_core_29_82_4 \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database solanum_tuberosum_core_29_82_4 does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Stub does not exist"
fi

