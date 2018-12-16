#!/usr/bin/env bash

# Arabidopsis lyrata
if [ -d /home/wangq/data/alignment/Ensembl/Alyr ]; then
    echo "==> Arabidopsis lyrata"

    if [ -f /home/wangq/data/alignment/Ensembl/Alyr/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Alyr/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Alyr

        find /home/wangq/data/ensembl94/gff3/arabidopsis_lyrata/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Alyr does not exist"
fi

# Arabidopsis thaliana
if [ -d /home/wangq/data/alignment/Ensembl/Atha ]; then
    echo "==> Arabidopsis thaliana"

    if [ -f /home/wangq/data/alignment/Ensembl/Atha/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Atha/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Atha

        find /home/wangq/data/ensembl94/gff3/arabidopsis_thaliana/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Atha does not exist"
fi

# Aspergillus fumigatus
if [ -d /home/wangq/data/alignment/Ensembl/Afum ]; then
    echo "==> Aspergillus fumigatus"

    if [ -f /home/wangq/data/alignment/Ensembl/Afum/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Afum/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Afum

        find /home/wangq/data/ensembl94/gff3/aspergillus_fumigatus/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Afum does not exist"
fi

# Brachypodium distachyon
if [ -d /home/wangq/data/alignment/Ensembl/Bdis ]; then
    echo "==> Brachypodium distachyon"

    if [ -f /home/wangq/data/alignment/Ensembl/Bdis/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Bdis/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Bdis

        find /home/wangq/data/ensembl94/gff3/brachypodium_distachyon/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Bdis does not exist"
fi

# Brassica oleracea
if [ -d /home/wangq/data/alignment/Ensembl/Bole ]; then
    echo "==> Brassica oleracea"

    if [ -f /home/wangq/data/alignment/Ensembl/Bole/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Bole/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Bole

        find /home/wangq/data/ensembl94/gff3/brassica_oleracea/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Bole does not exist"
fi

# Brassica rapa
if [ -d /home/wangq/data/alignment/Ensembl/Brap ]; then
    echo "==> Brassica rapa"

    if [ -f /home/wangq/data/alignment/Ensembl/Brap/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Brap/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Brap

        find /home/wangq/data/ensembl94/gff3/brassica_rapa/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Brap does not exist"
fi

# Caenorhabditis briggsae
if [ -d /home/wangq/data/alignment/Ensembl/Cbri ]; then
    echo "==> Caenorhabditis briggsae"

    if [ -f /home/wangq/data/alignment/Ensembl/Cbri/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Cbri/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Cbri

        find /home/wangq/data/ensembl94/gff3/caenorhabditis_briggsae/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Cbri does not exist"
fi

# Caenorhabditis elegans
if [ -d /home/wangq/data/alignment/Ensembl/Cele ]; then
    echo "==> Caenorhabditis elegans"

    if [ -f /home/wangq/data/alignment/Ensembl/Cele/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Cele/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Cele

        find /home/wangq/data/ensembl94/gff3/caenorhabditis_elegans/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Cele does not exist"
fi

# Dictyostelium discoideum
if [ -d /home/wangq/data/alignment/Ensembl/Ddis ]; then
    echo "==> Dictyostelium discoideum"

    if [ -f /home/wangq/data/alignment/Ensembl/Ddis/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Ddis/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Ddis

        find /home/wangq/data/ensembl94/gff3/dictyostelium_discoideum/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Ddis does not exist"
fi

# Drosophila melanogaster
if [ -d /home/wangq/data/alignment/Ensembl/Dmel ]; then
    echo "==> Drosophila melanogaster"

    if [ -f /home/wangq/data/alignment/Ensembl/Dmel/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Dmel/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Dmel

        find /home/wangq/data/ensembl94/gff3/drosophila_melanogaster/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Dmel does not exist"
fi

# Drosophila simulans
if [ -d /home/wangq/data/alignment/Ensembl/Dsim ]; then
    echo "==> Drosophila simulans"

    if [ -f /home/wangq/data/alignment/Ensembl/Dsim/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Dsim/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Dsim

        find /home/wangq/data/ensembl94/gff3/drosophila_simulans/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Dsim does not exist"
fi

# Glycine max
if [ -d /home/wangq/data/alignment/Ensembl/Gmax ]; then
    echo "==> Glycine max"

    if [ -f /home/wangq/data/alignment/Ensembl/Gmax/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Gmax/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Gmax

        find /home/wangq/data/ensembl94/gff3/glycine_max/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Gmax does not exist"
fi

# Gorilla gorilla
if [ -d /home/wangq/data/alignment/Ensembl/Gorilla ]; then
    echo "==> Gorilla gorilla"

    if [ -f /home/wangq/data/alignment/Ensembl/Gorilla/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Gorilla/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Gorilla

        find /home/wangq/data/ensembl94/gff3/gorilla_gorilla/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Gorilla does not exist"
fi

# Homo sapiens
if [ -d /home/wangq/data/alignment/Ensembl/Human ]; then
    echo "==> Homo sapiens"

    if [ -f /home/wangq/data/alignment/Ensembl/Human/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Human/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Human

        find /home/wangq/data/ensembl94/gff3/homo_sapiens/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Human does not exist"
fi

# Macaca mulatta
if [ -d /home/wangq/data/alignment/Ensembl/Rhesus ]; then
    echo "==> Macaca mulatta"

    if [ -f /home/wangq/data/alignment/Ensembl/Rhesus/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Rhesus/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Rhesus

        find /home/wangq/data/ensembl94/gff3/macaca_mulatta/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Rhesus does not exist"
fi

# Medicago truncatula
if [ -d /home/wangq/data/alignment/Ensembl/Mtru ]; then
    echo "==> Medicago truncatula"

    if [ -f /home/wangq/data/alignment/Ensembl/Mtru/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Mtru/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Mtru

        find /home/wangq/data/ensembl94/gff3/medicago_truncatula/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Mtru does not exist"
fi

# Mus musculus
if [ -d /home/wangq/data/alignment/Ensembl/Mouse ]; then
    echo "==> Mus musculus"

    if [ -f /home/wangq/data/alignment/Ensembl/Mouse/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Mouse/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Mouse

        find /home/wangq/data/ensembl94/gff3/mus_musculus/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Mouse does not exist"
fi

# Musa acuminata
if [ -d /home/wangq/data/alignment/Ensembl/Macu ]; then
    echo "==> Musa acuminata"

    if [ -f /home/wangq/data/alignment/Ensembl/Macu/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Macu/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Macu

        find /home/wangq/data/ensembl94/gff3/musa_acuminata/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Macu does not exist"
fi

# Oryza indica
if [ -d /home/wangq/data/alignment/Ensembl/OsatInd ]; then
    echo "==> Oryza indica"

    if [ -f /home/wangq/data/alignment/Ensembl/OsatInd/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/OsatInd/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/OsatInd

        find /home/wangq/data/ensembl94/gff3/oryza_indica/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatInd does not exist"
fi

# Oryza sativa
if [ -d /home/wangq/data/alignment/Ensembl/OsatJap ]; then
    echo "==> Oryza sativa"

    if [ -f /home/wangq/data/alignment/Ensembl/OsatJap/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/OsatJap/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/OsatJap

        find /home/wangq/data/ensembl94/gff3/oryza_sativa/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatJap does not exist"
fi

# Pan troglodytes
if [ -d /home/wangq/data/alignment/Ensembl/Chimp ]; then
    echo "==> Pan troglodytes"

    if [ -f /home/wangq/data/alignment/Ensembl/Chimp/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Chimp/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Chimp

        find /home/wangq/data/ensembl94/gff3/pan_troglodytes/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Chimp does not exist"
fi

# Plasmodium falciparum
if [ -d /home/wangq/data/alignment/Ensembl/Pfal ]; then
    echo "==> Plasmodium falciparum"

    if [ -f /home/wangq/data/alignment/Ensembl/Pfal/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Pfal/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Pfal

        find /home/wangq/data/ensembl94/gff3/plasmodium_falciparum/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Pfal does not exist"
fi

# Pongo abelii
if [ -d /home/wangq/data/alignment/Ensembl/Orangutan ]; then
    echo "==> Pongo abelii"

    if [ -f /home/wangq/data/alignment/Ensembl/Orangutan/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Orangutan/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Orangutan

        find /home/wangq/data/ensembl94/gff3/pongo_abelii/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Orangutan does not exist"
fi

# Rattus norvegicus
if [ -d /home/wangq/data/alignment/Ensembl/Rat ]; then
    echo "==> Rattus norvegicus"

    if [ -f /home/wangq/data/alignment/Ensembl/Rat/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Rat/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Rat

        find /home/wangq/data/ensembl94/gff3/rattus_norvegicus/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Rat does not exist"
fi

# Saccharomyces cerevisiae
if [ -d /home/wangq/data/alignment/Ensembl/S288c ]; then
    echo "==> Saccharomyces cerevisiae"

    if [ -f /home/wangq/data/alignment/Ensembl/S288c/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/S288c/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/S288c

        find /home/wangq/data/ensembl94/gff3/saccharomyces_cerevisiae/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/S288c does not exist"
fi

# Schizosaccharomyces pombe
if [ -d /home/wangq/data/alignment/Ensembl/Spom ]; then
    echo "==> Schizosaccharomyces pombe"

    if [ -f /home/wangq/data/alignment/Ensembl/Spom/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Spom/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Spom

        find /home/wangq/data/ensembl94/gff3/schizosaccharomyces_pombe/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Spom does not exist"
fi

# Setaria italica
if [ -d /home/wangq/data/alignment/Ensembl/Sita ]; then
    echo "==> Setaria italica"

    if [ -f /home/wangq/data/alignment/Ensembl/Sita/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Sita/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Sita

        find /home/wangq/data/ensembl94/gff3/setaria_italica/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Sita does not exist"
fi

# Solanum lycopersicum
if [ -d /home/wangq/data/alignment/Ensembl/Slyc ]; then
    echo "==> Solanum lycopersicum"

    if [ -f /home/wangq/data/alignment/Ensembl/Slyc/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Slyc/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Slyc

        find /home/wangq/data/ensembl94/gff3/solanum_lycopersicum/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Slyc does not exist"
fi

# Solanum tuberosum
if [ -d /home/wangq/data/alignment/Ensembl/Stub ]; then
    echo "==> Solanum tuberosum"

    if [ -f /home/wangq/data/alignment/Ensembl/Stub/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Stub/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Stub

        find /home/wangq/data/ensembl94/gff3/solanum_tuberosum/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Stub does not exist"
fi

# Sorghum bicolor
if [ -d /home/wangq/data/alignment/Ensembl/Sbic ]; then
    echo "==> Sorghum bicolor"

    if [ -f /home/wangq/data/alignment/Ensembl/Sbic/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Sbic/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Sbic

        find /home/wangq/data/ensembl94/gff3/sorghum_bicolor/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Sbic does not exist"
fi

# Vitis vinifera
if [ -d /home/wangq/data/alignment/Ensembl/Vvin ]; then
    echo "==> Vitis vinifera"

    if [ -f /home/wangq/data/alignment/Ensembl/Vvin/anno.yml ]; then
        echo "==> /home/wangq/data/alignment/Ensembl/Vvin/anno.yml exists"
    else
        cd /home/wangq/data/alignment/Ensembl/Vvin

        find /home/wangq/data/ensembl94/gff3/vitis_vinifera/ -name "*.gff3.gz" |
            grep -v "abinitio.gff3" |
            grep -v "chr.gff3" |
            xargs gzip -d -c > chr.gff
        runlist gff --tag CDS --remove chr.gff -o cds.yml

        faops masked *.fa |
            jrunlist cover stdin -o repeat.yml

        runlist merge repeat.yml cds.yml -o anno.yml
        rm repeat.yml cds.yml
    fi
else
    echo "==> /home/wangq/data/alignment/Ensembl/Vvin does not exist"
fi

