#/bin/bash

# Arabidopsis lyrata
if [ ! -d /home/wangq/data/alignment/Ensembl/Alyr ]
then
    echo "==> Arabidopsis lyrata"

    mkdir -p /home/wangq/data/alignment/Ensembl/Alyr
    cd /home/wangq/data/alignment/Ensembl/Alyr
    
    find /home/wangq/data/ensembl82/fasta/arabidopsis_lyrata/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat scaffold*.fa > Un.fa
    rm scaffold*.fa
    mv Un.fa Un.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/Alyr exists"
fi

# Arabidopsis thaliana
if [ ! -d /home/wangq/data/alignment/Ensembl/Atha ]
then
    echo "==> Arabidopsis thaliana"

    mkdir -p /home/wangq/data/alignment/Ensembl/Atha
    cd /home/wangq/data/alignment/Ensembl/Atha
    
    find /home/wangq/data/ensembl82/fasta/arabidopsis_thaliana/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    mv Mt.fa Mt.fa.skip
    mv Pt.fa Pt.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/Atha exists"
fi

# Brachypodium distachyon
if [ ! -d /home/wangq/data/alignment/Ensembl/Bdis ]
then
    echo "==> Brachypodium distachyon"

    mkdir -p /home/wangq/data/alignment/Ensembl/Bdis
    cd /home/wangq/data/alignment/Ensembl/Bdis
    
    find /home/wangq/data/ensembl82/fasta/brachypodium_distachyon/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat GG*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm GG*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Bdis exists"
fi

# Brassica oleracea
if [ ! -d /home/wangq/data/alignment/Ensembl/Bole ]
then
    echo "==> Brassica oleracea"

    mkdir -p /home/wangq/data/alignment/Ensembl/Bole
    cd /home/wangq/data/alignment/Ensembl/Bole
    
    find /home/wangq/data/ensembl82/fasta/brassica_oleracea/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat Scaffold*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm Scaffold*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Bole exists"
fi

# Brassica rapa
if [ ! -d /home/wangq/data/alignment/Ensembl/Brap ]
then
    echo "==> Brassica rapa"

    mkdir -p /home/wangq/data/alignment/Ensembl/Brap
    cd /home/wangq/data/alignment/Ensembl/Brap
    
    find /home/wangq/data/ensembl82/fasta/brassica_rapa/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat Scaffold*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm Scaffold*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Brap exists"
fi

# Caenorhabditis elegans
if [ ! -d /home/wangq/data/alignment/Ensembl/Cele ]
then
    echo "==> Caenorhabditis elegans"

    mkdir -p /home/wangq/data/alignment/Ensembl/Cele
    cd /home/wangq/data/alignment/Ensembl/Cele
    
    find /home/wangq/data/ensembl82/fasta/caenorhabditis_elegans/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    mv MtDNA.fa MtDNA.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/Cele exists"
fi

# Drosophila melanogaster
if [ ! -d /home/wangq/data/alignment/Ensembl/Dmel ]
then
    echo "==> Drosophila melanogaster"

    mkdir -p /home/wangq/data/alignment/Ensembl/Dmel
    cd /home/wangq/data/alignment/Ensembl/Dmel
    
    find /home/wangq/data/ensembl82/fasta/drosophila_melanogaster/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    rm *Scaffold*.fa 211*.fa
    mv 4.fa 4.fa.skip
    mv Y.fa Y.fa.skip
    mv rDNA.fa rDNA.fa.skip
    mv dmel_mitochondrion_genome.fa dmel_mitochondrion_genome.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/Dmel exists"
fi

# Drosophila simulans
if [ ! -d /home/wangq/data/alignment/Ensembl/Dsim ]
then
    echo "==> Drosophila simulans"

    mkdir -p /home/wangq/data/alignment/Ensembl/Dsim
    cd /home/wangq/data/alignment/Ensembl/Dsim
    
    find /home/wangq/data/ensembl82/fasta/drosophila_simulans/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile


else
    echo "==> /home/wangq/data/alignment/Ensembl/Dsim exists"
fi

# Gorilla gorilla
if [ ! -d /home/wangq/data/alignment/Ensembl/Gorilla ]
then
    echo "==> Gorilla gorilla"

    mkdir -p /home/wangq/data/alignment/Ensembl/Gorilla
    cd /home/wangq/data/alignment/Ensembl/Gorilla
    
    find /home/wangq/data/ensembl82/fasta/gorilla_gorilla/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile


else
    echo "==> /home/wangq/data/alignment/Ensembl/Gorilla exists"
fi

# Homo sapiens
if [ ! -d /home/wangq/data/alignment/Ensembl/Human ]
then
    echo "==> Homo sapiens"

    mkdir -p /home/wangq/data/alignment/Ensembl/Human
    cd /home/wangq/data/alignment/Ensembl/Human
    
    find /home/wangq/data/ensembl82/fasta/homo_sapiens/dna/ -name "*dna_sm.primary_assembly*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat GL*.fa > Un.fa
    rm GL*.fa
    mv Un.fa Un.fa.skip
    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/Human exists"
fi

# Macaca mulatta
if [ ! -d /home/wangq/data/alignment/Ensembl/Rhesus ]
then
    echo "==> Macaca mulatta"

    mkdir -p /home/wangq/data/alignment/Ensembl/Rhesus
    cd /home/wangq/data/alignment/Ensembl/Rhesus
    
    find /home/wangq/data/ensembl82/fasta/macaca_mulatta/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat 1099*.fa > Un.fa
    rm 1099*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Rhesus exists"
fi

# Mus musculus
if [ ! -d /home/wangq/data/alignment/Ensembl/Mouse ]
then
    echo "==> Mus musculus"

    mkdir -p /home/wangq/data/alignment/Ensembl/Mouse
    cd /home/wangq/data/alignment/Ensembl/Mouse
    
    find /home/wangq/data/ensembl82/fasta/mus_musculus/dna/ -name "*dna_sm.primary_assembly*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

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
if [ ! -d /home/wangq/data/alignment/Ensembl/Macu ]
then
    echo "==> Musa acuminata"

    mkdir -p /home/wangq/data/alignment/Ensembl/Macu
    cd /home/wangq/data/alignment/Ensembl/Macu
    
    find /home/wangq/data/ensembl82/fasta/musa_acuminata/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile


else
    echo "==> /home/wangq/data/alignment/Ensembl/Macu exists"
fi

# Oryza indica
if [ ! -d /home/wangq/data/alignment/Ensembl/OsatInd ]
then
    echo "==> Oryza indica"

    mkdir -p /home/wangq/data/alignment/Ensembl/OsatInd
    cd /home/wangq/data/alignment/Ensembl/OsatInd
    
    find /home/wangq/data/ensembl82/fasta/oryza_indica/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat AA*.fa > Un.fa
    cat CH*.fa >> Un.fa
    cat Sup*.fa >> Un.fa
    rm AA*.fa CH*.fa Sup*.fa
    mv Un.fa Un.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatInd exists"
fi

# Oryza sativa
if [ ! -d /home/wangq/data/alignment/Ensembl/OsatJap ]
then
    echo "==> Oryza sativa"

    mkdir -p /home/wangq/data/alignment/Ensembl/OsatJap
    cd /home/wangq/data/alignment/Ensembl/OsatJap
    
    find /home/wangq/data/ensembl82/fasta/oryza_sativa/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat AC*.fa > Un.fa
    cat AP*.fa >> Un.fa
    cat Syng*.fa >> Un.fa
    rm AC*.fa AP*.fa Syng*.fa
    mv Un.fa Un.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/OsatJap exists"
fi

# Pan troglodytes
if [ ! -d /home/wangq/data/alignment/Ensembl/Chimp ]
then
    echo "==> Pan troglodytes"

    mkdir -p /home/wangq/data/alignment/Ensembl/Chimp
    cd /home/wangq/data/alignment/Ensembl/Chimp
    
    find /home/wangq/data/ensembl82/fasta/pan_troglodytes/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat GL*.fa > Un.fa
    cat AACZ*.fa >> Un.fa
    rm GL*.fa AACZ*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Chimp exists"
fi

# Pongo abelii
if [ ! -d /home/wangq/data/alignment/Ensembl/Orangutan ]
then
    echo "==> Pongo abelii"

    mkdir -p /home/wangq/data/alignment/Ensembl/Orangutan
    cd /home/wangq/data/alignment/Ensembl/Orangutan
    
    find /home/wangq/data/ensembl82/fasta/pongo_abelii/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile


else
    echo "==> /home/wangq/data/alignment/Ensembl/Orangutan exists"
fi

# Rattus norvegicus
if [ ! -d /home/wangq/data/alignment/Ensembl/Rat ]
then
    echo "==> Rattus norvegicus"

    mkdir -p /home/wangq/data/alignment/Ensembl/Rat
    cd /home/wangq/data/alignment/Ensembl/Rat
    
    find /home/wangq/data/ensembl82/fasta/rattus_norvegicus/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

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
if [ ! -d /home/wangq/data/alignment/Ensembl/S288c ]
then
    echo "==> Saccharomyces cerevisiae"

    mkdir -p /home/wangq/data/alignment/Ensembl/S288c
    cd /home/wangq/data/alignment/Ensembl/S288c
    
    find /home/wangq/data/ensembl82/fasta/saccharomyces_cerevisiae/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    mv Mito.fa Mito.fa.skip
else
    echo "==> /home/wangq/data/alignment/Ensembl/S288c exists"
fi

# Setaria italica
if [ ! -d /home/wangq/data/alignment/Ensembl/Sita ]
then
    echo "==> Setaria italica"

    mkdir -p /home/wangq/data/alignment/Ensembl/Sita
    cd /home/wangq/data/alignment/Ensembl/Sita
    
    find /home/wangq/data/ensembl82/fasta/setaria_italica/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat Scaffold*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm Scaffold*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Sita exists"
fi

# Solanum lycopersicum
if [ ! -d /home/wangq/data/alignment/Ensembl/Slyc ]
then
    echo "==> Solanum lycopersicum"

    mkdir -p /home/wangq/data/alignment/Ensembl/Slyc
    cd /home/wangq/data/alignment/Ensembl/Slyc
    
    find /home/wangq/data/ensembl82/fasta/solanum_lycopersicum/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat SL2*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm SL2*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Slyc exists"
fi

# Solanum tuberosum
if [ ! -d /home/wangq/data/alignment/Ensembl/Stub ]
then
    echo "==> Solanum tuberosum"

    mkdir -p /home/wangq/data/alignment/Ensembl/Stub
    cd /home/wangq/data/alignment/Ensembl/Stub
    
    find /home/wangq/data/ensembl82/fasta/solanum_tuberosum/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile


else
    echo "==> /home/wangq/data/alignment/Ensembl/Stub exists"
fi

# Sorghum bicolor
if [ ! -d /home/wangq/data/alignment/Ensembl/Sbic ]
then
    echo "==> Sorghum bicolor"

    mkdir -p /home/wangq/data/alignment/Ensembl/Sbic
    cd /home/wangq/data/alignment/Ensembl/Sbic
    
    find /home/wangq/data/ensembl82/fasta/sorghum_bicolor/dna/ -name "*dna_sm.toplevel*" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat GL*.fa > Un.fa
    mv Un.fa Un.fa.skip
    rm GL*.fa
else
    echo "==> /home/wangq/data/alignment/Ensembl/Sbic exists"
fi

