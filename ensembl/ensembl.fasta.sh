#!/usr/bin/env bash

# Arabidopsis lyrata
if [ ! -d /home/wangq/data/alignment/Ensembl/Alyr ]; then
    echo "==> Arabidopsis lyrata"

    mkdir -p /home/wangq/data/alignment/Ensembl/Alyr
    cd /home/wangq/data/alignment/Ensembl/Alyr
    
    find /home/wangq/data/ensembl94/fasta/arabidopsis_lyrata/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat scaffold*.fa > Un.fa
    rm scaffold*.fa
    mv Un.fa Un.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Alyr exists"
fi

# Arabidopsis thaliana
if [ ! -d /home/wangq/data/alignment/Ensembl/Atha ]; then
    echo "==> Arabidopsis thaliana"

    mkdir -p /home/wangq/data/alignment/Ensembl/Atha
    cd /home/wangq/data/alignment/Ensembl/Atha
    
    find /home/wangq/data/ensembl94/fasta/arabidopsis_thaliana/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    mv Mt.fa Mt.fa.skip
    mv Pt.fa Pt.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Atha exists"
fi

# Aspergillus fumigatus
if [ ! -d /home/wangq/data/alignment/Ensembl/Afum ]; then
    echo "==> Aspergillus fumigatus"

    mkdir -p /home/wangq/data/alignment/Ensembl/Afum
    cd /home/wangq/data/alignment/Ensembl/Afum
    
    find /home/wangq/data/ensembl94/fasta/aspergillus_fumigatus/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile



else
    echo "==> /home/wangq/data/alignment/Ensembl/Afum exists"
fi

# Brachypodium distachyon
if [ ! -d /home/wangq/data/alignment/Ensembl/Bdis ]; then
    echo "==> Brachypodium distachyon"

    mkdir -p /home/wangq/data/alignment/Ensembl/Bdis
    cd /home/wangq/data/alignment/Ensembl/Bdis
    
    find /home/wangq/data/ensembl94/fasta/brachypodium_distachyon/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat KZ*.fa Bd1_*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm KZ*.fa Bd1_*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Bdis exists"
fi

# Brassica oleracea
if [ ! -d /home/wangq/data/alignment/Ensembl/Bole ]; then
    echo "==> Brassica oleracea"

    mkdir -p /home/wangq/data/alignment/Ensembl/Bole
    cd /home/wangq/data/alignment/Ensembl/Bole
    
    find /home/wangq/data/ensembl94/fasta/brassica_oleracea/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat Scaffold*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm Scaffold*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Bole exists"
fi

# Brassica rapa
if [ ! -d /home/wangq/data/alignment/Ensembl/Brap ]; then
    echo "==> Brassica rapa"

    mkdir -p /home/wangq/data/alignment/Ensembl/Brap
    cd /home/wangq/data/alignment/Ensembl/Brap
    
    find /home/wangq/data/ensembl94/fasta/brassica_rapa/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat Scaffold*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm Scaffold*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Brap exists"
fi

# Caenorhabditis briggsae
if [ ! -d /home/wangq/data/alignment/Ensembl/Cbri ]; then
    echo "==> Caenorhabditis briggsae"

    mkdir -p /home/wangq/data/alignment/Ensembl/Cbri
    cd /home/wangq/data/alignment/Ensembl/Cbri
    
    find /home/wangq/data/ensembl94/fasta/caenorhabditis_briggsae/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat cb25*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm cb25*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Cbri exists"
fi

# Caenorhabditis elegans
if [ ! -d /home/wangq/data/alignment/Ensembl/Cele ]; then
    echo "==> Caenorhabditis elegans"

    mkdir -p /home/wangq/data/alignment/Ensembl/Cele
    cd /home/wangq/data/alignment/Ensembl/Cele
    
    find /home/wangq/data/ensembl94/fasta/caenorhabditis_elegans/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    mv MtDNA.fa MtDNA.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Cele exists"
fi

# Dictyostelium discoideum
if [ ! -d /home/wangq/data/alignment/Ensembl/Ddis ]; then
    echo "==> Dictyostelium discoideum"

    mkdir -p /home/wangq/data/alignment/Ensembl/Ddis
    cd /home/wangq/data/alignment/Ensembl/Ddis
    
    find /home/wangq/data/ensembl94/fasta/dictyostelium_discoideum/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat CH*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm CH*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Ddis exists"
fi

# Drosophila melanogaster
if [ ! -d /home/wangq/data/alignment/Ensembl/Dmel ]; then
    echo "==> Drosophila melanogaster"

    mkdir -p /home/wangq/data/alignment/Ensembl/Dmel
    cd /home/wangq/data/alignment/Ensembl/Dmel
    
    find /home/wangq/data/ensembl94/fasta/drosophila_melanogaster/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat *Scaffold*.fa 211*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm *Scaffold*.fa 211*.fa
    mv 4.fa 4.fa.skip
    mv Y.fa Y.fa.skip
    mv rDNA.fa rDNA.fa.skip
    mv mitochondrion_genome.fa mitochondrion_genome.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Dmel exists"
fi

# Drosophila simulans
if [ ! -d /home/wangq/data/alignment/Ensembl/Dsim ]; then
    echo "==> Drosophila simulans"

    mkdir -p /home/wangq/data/alignment/Ensembl/Dsim
    cd /home/wangq/data/alignment/Ensembl/Dsim
    
    find /home/wangq/data/ensembl94/fasta/drosophila_simulans/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat JPYS*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm JPYS*.fa
    mv 4.fa 4.fa.skip
    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Dsim exists"
fi

# Glycine max
if [ ! -d /home/wangq/data/alignment/Ensembl/Gmax ]; then
    echo "==> Glycine max"

    mkdir -p /home/wangq/data/alignment/Ensembl/Gmax
    cd /home/wangq/data/alignment/Ensembl/Gmax
    
    find /home/wangq/data/ensembl94/fasta/glycine_max/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat scaffold*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm scaffold*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Gmax exists"
fi

# Gorilla gorilla
if [ ! -d /home/wangq/data/alignment/Ensembl/Gorilla ]; then
    echo "==> Gorilla gorilla"

    mkdir -p /home/wangq/data/alignment/Ensembl/Gorilla
    cd /home/wangq/data/alignment/Ensembl/Gorilla
    
    find /home/wangq/data/ensembl94/fasta/gorilla_gorilla/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat CABD*.fa > Un.fa
    rm CABD*.fa
    mv Un.fa Un.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Gorilla exists"
fi

# Homo sapiens
if [ ! -d /home/wangq/data/alignment/Ensembl/Human ]; then
    echo "==> Homo sapiens"

    mkdir -p /home/wangq/data/alignment/Ensembl/Human
    cd /home/wangq/data/alignment/Ensembl/Human
    
    find /home/wangq/data/ensembl94/fasta/homo_sapiens/dna/ -name "*dna_sm.primary_assembly*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat GL*.fa > Un.fa
    rm GL*.fa
    mv Un.fa Un.fa.skip
    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Human exists"
fi

# Macaca mulatta
if [ ! -d /home/wangq/data/alignment/Ensembl/Rhesus ]; then
    echo "==> Macaca mulatta"

    mkdir -p /home/wangq/data/alignment/Ensembl/Rhesus
    cd /home/wangq/data/alignment/Ensembl/Rhesus
    
    find /home/wangq/data/ensembl94/fasta/macaca_mulatta/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat JSUE*.fa KQ*.fa > Un.fa
    rm JSUE*.fa KQ*.fa
    mv Un.fa Un.fa.skip
    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Rhesus exists"
fi

# Medicago truncatula
if [ ! -d /home/wangq/data/alignment/Ensembl/Mtru ]; then
    echo "==> Medicago truncatula"

    mkdir -p /home/wangq/data/alignment/Ensembl/Mtru
    cd /home/wangq/data/alignment/Ensembl/Mtru
    
    find /home/wangq/data/ensembl94/fasta/medicago_truncatula/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat scaffold*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm scaffold*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Mtru exists"
fi

# Mus musculus
if [ ! -d /home/wangq/data/alignment/Ensembl/Mouse ]; then
    echo "==> Mus musculus"

    mkdir -p /home/wangq/data/alignment/Ensembl/Mouse
    cd /home/wangq/data/alignment/Ensembl/Mouse
    
    find /home/wangq/data/ensembl94/fasta/mus_musculus/dna/ -name "*dna_sm.primary_assembly*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat GL*.fa > Un.fa
    cat JH*.fa >> Un.fa
    rm GL*.fa JH*.fa
    mv Un.fa Un.fa.skip
    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Mouse exists"
fi

# Musa acuminata
if [ ! -d /home/wangq/data/alignment/Ensembl/Macu ]; then
    echo "==> Musa acuminata"

    mkdir -p /home/wangq/data/alignment/Ensembl/Macu
    cd /home/wangq/data/alignment/Ensembl/Macu
    
    find /home/wangq/data/ensembl94/fasta/musa_acuminata/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile



else
    echo "==> /home/wangq/data/alignment/Ensembl/Macu exists"
fi

# Oryza indica
if [ ! -d /home/wangq/data/alignment/Ensembl/OsatInd ]; then
    echo "==> Oryza indica"

    mkdir -p /home/wangq/data/alignment/Ensembl/OsatInd
    cd /home/wangq/data/alignment/Ensembl/OsatInd
    
    find /home/wangq/data/ensembl94/fasta/oryza_indica/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat AA*.fa > Un.fa
    cat CH*.fa >> Un.fa
    cat Sup*.fa >> Un.fa
    rm AA*.fa CH*.fa Sup*.fa
    mv Un.fa Un.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatInd exists"
fi

# Oryza sativa
if [ ! -d /home/wangq/data/alignment/Ensembl/OsatJap ]; then
    echo "==> Oryza sativa"

    mkdir -p /home/wangq/data/alignment/Ensembl/OsatJap
    cd /home/wangq/data/alignment/Ensembl/OsatJap
    
    find /home/wangq/data/ensembl94/fasta/oryza_sativa/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat AC*.fa > Un.fa
    cat AP*.fa >> Un.fa
    cat Syng*.fa >> Un.fa
    rm AC*.fa AP*.fa Syng*.fa
    mv Un.fa Un.fa.skip
    mv Mt.fa Mt.fa.skip
    mv Pt.fa Pt.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatJap exists"
fi

# Pan troglodytes
if [ ! -d /home/wangq/data/alignment/Ensembl/Chimp ]; then
    echo "==> Pan troglodytes"

    mkdir -p /home/wangq/data/alignment/Ensembl/Chimp
    cd /home/wangq/data/alignment/Ensembl/Chimp
    
    find /home/wangq/data/ensembl94/fasta/pan_troglodytes/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat KV*.fa AACZ*.fa > Un.fa
    rm KV*.fa AACZ*.fa
    mv Un.fa Un.fa.skip
    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Chimp exists"
fi

# Plasmodium falciparum
if [ ! -d /home/wangq/data/alignment/Ensembl/Pfal ]; then
    echo "==> Plasmodium falciparum"

    mkdir -p /home/wangq/data/alignment/Ensembl/Pfal
    cd /home/wangq/data/alignment/Ensembl/Pfal
    
    find /home/wangq/data/ensembl94/fasta/plasmodium_falciparum/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile



else
    echo "==> /home/wangq/data/alignment/Ensembl/Pfal exists"
fi

# Pongo abelii
if [ ! -d /home/wangq/data/alignment/Ensembl/Orangutan ]; then
    echo "==> Pongo abelii"

    mkdir -p /home/wangq/data/alignment/Ensembl/Orangutan
    cd /home/wangq/data/alignment/Ensembl/Orangutan
    
    find /home/wangq/data/ensembl94/fasta/pongo_abelii/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat *random.fa *hap*.fa >> Un.fa
    rm *random.fa *hap*.fa
    mv Un.fa Un.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Orangutan exists"
fi

# Rattus norvegicus
if [ ! -d /home/wangq/data/alignment/Ensembl/Rat ]; then
    echo "==> Rattus norvegicus"

    mkdir -p /home/wangq/data/alignment/Ensembl/Rat
    cd /home/wangq/data/alignment/Ensembl/Rat
    
    find /home/wangq/data/ensembl94/fasta/rattus_norvegicus/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat KL*.fa > Un.fa
    cat AABR*.fa >> Un.fa
    rm KL*.fa AABR*.fa
    mv Un.fa Un.fa.skip
    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Rat exists"
fi

# Saccharomyces cerevisiae
if [ ! -d /home/wangq/data/alignment/Ensembl/S288c ]; then
    echo "==> Saccharomyces cerevisiae"

    mkdir -p /home/wangq/data/alignment/Ensembl/S288c
    cd /home/wangq/data/alignment/Ensembl/S288c
    
    find /home/wangq/data/ensembl94/fasta/saccharomyces_cerevisiae/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    mv Mito.fa Mito.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/S288c exists"
fi

# Schizosaccharomyces pombe
if [ ! -d /home/wangq/data/alignment/Ensembl/Spom ]; then
    echo "==> Schizosaccharomyces pombe"

    mkdir -p /home/wangq/data/alignment/Ensembl/Spom
    cd /home/wangq/data/alignment/Ensembl/Spom
    
    find /home/wangq/data/ensembl94/fasta/schizosaccharomyces_pombe/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile



else
    echo "==> /home/wangq/data/alignment/Ensembl/Spom exists"
fi

# Setaria italica
if [ ! -d /home/wangq/data/alignment/Ensembl/Sita ]; then
    echo "==> Setaria italica"

    mkdir -p /home/wangq/data/alignment/Ensembl/Sita
    cd /home/wangq/data/alignment/Ensembl/Sita
    
    find /home/wangq/data/ensembl94/fasta/setaria_italica/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat KQ*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm KQ*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Sita exists"
fi

# Solanum lycopersicum
if [ ! -d /home/wangq/data/alignment/Ensembl/Slyc ]; then
    echo "==> Solanum lycopersicum"

    mkdir -p /home/wangq/data/alignment/Ensembl/Slyc
    cd /home/wangq/data/alignment/Ensembl/Slyc
    
    find /home/wangq/data/ensembl94/fasta/solanum_lycopersicum/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat SL2*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm SL2*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Slyc exists"
fi

# Solanum tuberosum
if [ ! -d /home/wangq/data/alignment/Ensembl/Stub ]; then
    echo "==> Solanum tuberosum"

    mkdir -p /home/wangq/data/alignment/Ensembl/Stub
    cd /home/wangq/data/alignment/Ensembl/Stub
    
    find /home/wangq/data/ensembl94/fasta/solanum_tuberosum/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile



else
    echo "==> /home/wangq/data/alignment/Ensembl/Stub exists"
fi

# Sorghum bicolor
if [ ! -d /home/wangq/data/alignment/Ensembl/Sbic ]; then
    echo "==> Sorghum bicolor"

    mkdir -p /home/wangq/data/alignment/Ensembl/Sbic
    cd /home/wangq/data/alignment/Ensembl/Sbic
    
    find /home/wangq/data/ensembl94/fasta/sorghum_bicolor/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat super*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm super*.fa

else
    echo "==> /home/wangq/data/alignment/Ensembl/Sbic exists"
fi

# Vitis vinifera
if [ ! -d /home/wangq/data/alignment/Ensembl/Vvin ]; then
    echo "==> Vitis vinifera"

    mkdir -p /home/wangq/data/alignment/Ensembl/Vvin
    cd /home/wangq/data/alignment/Ensembl/Vvin
    
    find /home/wangq/data/ensembl94/fasta/vitis_vinifera/dna/ -name "*dna_sm.toplevel*" |
        xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa |
        perl -nla -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
        ' |
        uniq > listFile
    faops some toplevel.fa listFile stdout |
        faops filter -N stdin stdout |
        faops split-name stdin .
    rm toplevel.fa listFile

    cat *random.fa >> Un.fa
    rm *random.fa
    mv Un.fa Un.fa.skip

else
    echo "==> /home/wangq/data/alignment/Ensembl/Vvin exists"
fi

