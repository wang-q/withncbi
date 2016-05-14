#!/usr/bin/env bash

# Arabidopsis thaliana
if [ -d /home/wangq/data/alignment/Ensembl/Atha ];
then
    echo "==> Arabidopsis thaliana"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Arabidopsis_thaliana';
    then
        cd /home/wangq/data/alignment/Ensembl/Atha

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Arabidopsis_thaliana \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Arabidopsis_thaliana \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Arabidopsis_thaliana does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Atha does not exist"
fi

# Brassica rapa
if [ -d /home/wangq/data/alignment/Ensembl/Brap ];
then
    echo "==> Brassica rapa"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Brassica_rapa';
    then
        cd /home/wangq/data/alignment/Ensembl/Brap

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Brassica_rapa \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Brassica_rapa \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Brassica_rapa does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Brap does not exist"
fi

# Caenorhabditis elegans
if [ -d /home/wangq/data/alignment/Ensembl/Cele ];
then
    echo "==> Caenorhabditis elegans"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Caenorhabditis_elegans';
    then
        cd /home/wangq/data/alignment/Ensembl/Cele

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Caenorhabditis_elegans \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Caenorhabditis_elegans \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Caenorhabditis_elegans does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Cele does not exist"
fi

# Drosophila melanogaster
if [ -d /home/wangq/data/alignment/Ensembl/Dmel ];
then
    echo "==> Drosophila melanogaster"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Drosophila_melanogaster';
    then
        cd /home/wangq/data/alignment/Ensembl/Dmel

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Drosophila_melanogaster \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Drosophila_melanogaster \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Drosophila_melanogaster does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Dmel does not exist"
fi

# Drosophila simulans
if [ -d /home/wangq/data/alignment/Ensembl/Dsim ];
then
    echo "==> Drosophila simulans"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Drosophila_simulans';
    then
        cd /home/wangq/data/alignment/Ensembl/Dsim

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Drosophila_simulans \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Drosophila_simulans \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Drosophila_simulans does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Dsim does not exist"
fi

# Gorilla gorilla
if [ -d /home/wangq/data/alignment/Ensembl/Gorilla ];
then
    echo "==> Gorilla gorilla"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Gorilla_gorilla';
    then
        cd /home/wangq/data/alignment/Ensembl/Gorilla

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Gorilla_gorilla \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Gorilla_gorilla \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Gorilla_gorilla does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Gorilla does not exist"
fi

# Homo sapiens
if [ -d /home/wangq/data/alignment/Ensembl/Human ];
then
    echo "==> Homo sapiens"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Homo_sapiens';
    then
        cd /home/wangq/data/alignment/Ensembl/Human

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Homo_sapiens \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Homo_sapiens \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Homo_sapiens does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Human does not exist"
fi

# Macaca mulatta
if [ -d /home/wangq/data/alignment/Ensembl/Rhesus ];
then
    echo "==> Macaca mulatta"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Macaca_mulatta';
    then
        cd /home/wangq/data/alignment/Ensembl/Rhesus

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Macaca_mulatta \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Macaca_mulatta \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Macaca_mulatta does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Rhesus does not exist"
fi

# Mus musculus
if [ -d /home/wangq/data/alignment/Ensembl/Mouse ];
then
    echo "==> Mus musculus"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Mus_musculus';
    then
        cd /home/wangq/data/alignment/Ensembl/Mouse

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Mus_musculus \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Mus_musculus \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Mus_musculus does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Mouse does not exist"
fi

# Oryza sativa
if [ -d /home/wangq/data/alignment/Ensembl/OsatJap ];
then
    echo "==> Oryza sativa"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Oryza_sativa';
    then
        cd /home/wangq/data/alignment/Ensembl/OsatJap

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Oryza_sativa \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Oryza_sativa \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Oryza_sativa does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatJap does not exist"
fi

# Pan troglodytes
if [ -d /home/wangq/data/alignment/Ensembl/Chimp ];
then
    echo "==> Pan troglodytes"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Pan_troglodytes';
    then
        cd /home/wangq/data/alignment/Ensembl/Chimp

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Pan_troglodytes \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Pan_troglodytes \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Pan_troglodytes does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Chimp does not exist"
fi

# Pongo abelii
if [ -d /home/wangq/data/alignment/Ensembl/Orangutan ];
then
    echo "==> Pongo abelii"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Pongo_abelii';
    then
        cd /home/wangq/data/alignment/Ensembl/Orangutan

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Pongo_abelii \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Pongo_abelii \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Pongo_abelii does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Orangutan does not exist"
fi

# Saccharomyces cerevisiae
if [ -d /home/wangq/data/alignment/Ensembl/S288c ];
then
    echo "==> Saccharomyces cerevisiae"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Saccharomyces_cerevisiae';
    then
        cd /home/wangq/data/alignment/Ensembl/S288c

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Saccharomyces_cerevisiae \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Saccharomyces_cerevisiae \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Saccharomyces_cerevisiae does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/S288c does not exist"
fi

# Solanum lycopersicum
if [ -d /home/wangq/data/alignment/Ensembl/Slyc ];
then
    echo "==> Solanum lycopersicum"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Solanum_lycopersicum';
    then
        cd /home/wangq/data/alignment/Ensembl/Slyc

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Solanum_lycopersicum \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Solanum_lycopersicum \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Solanum_lycopersicum does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Slyc does not exist"
fi

# Solanum tuberosum
if [ -d /home/wangq/data/alignment/Ensembl/Stub ];
then
    echo "==> Solanum tuberosum"

    if mysql -hlocalhost -ualignDB -palignDB -e 'use Solanum_tuberosum';
    then
        cd /home/wangq/data/alignment/Ensembl/Stub

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Solanum_tuberosum \
            -f repeat \
            -o repeat.yml

        perl ~/Scripts/withncbi/ensembl/feature_runlists.pl \
            -e Solanum_tuberosum \
            -f cds \
            -o cds.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    else
        echo "==> Database Solanum_tuberosum does not exist"
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Stub does not exist"
fi

