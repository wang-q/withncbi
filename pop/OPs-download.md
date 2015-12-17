# Operating steps for each groups

Less detailed than Trichoderma in [README.md](README.md), but include examples
for genomes out of WGS, which usually in better assembling levels.

## *Saccharomyces* WGS

1. Create `pop/saccharomyces.tsv` manually.

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Saccharomyces&order=organism
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid4930[Organism:exp]

    ```bash
    export GENUS_ID=4930
    export GENUS=saccharomyces
    mkdir -p ~/data/alignment/Fungi/$GENUS          # operation directory
    mkdir -p ~/data/alignment/Fungi/GENOMES/$GENUS  # sequence directory

    cd ~/data/alignment/Fungi/GENOMES/$GENUS

    ... # paste codes from README.md

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

    Remove all strains of Saccharomyces cerevisiae. S288c will be injected later.

    ```bash
    mv saccharomyces.tsv all.tsv
    cat all.tsv | grep -v Scer_ > saccharomyces.tsv
    ```

    Edit them to fix names and comment out bad strains.

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/saccharomyces
    cd ~/data/alignment/Fungi/GENOMES/saccharomyces

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/saccharomyces.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/saccharomyces.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download *Saccharomyces cerevisiae* S288c.
  This step is totally manual operation. **Be careful.**

    | assigned name | organism_name                    | assembly_accession |
    | :------------ | :------------                    | :------------      |
    | Scer_S288c    | *Saccharomyces cerevisiae* S288c | GCF_000146045.2    |

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/saccharomyces/DOWNLOAD
    cd ~/data/alignment/Fungi/GENOMES/saccharomyces/DOWNLOAD

    # Omit chrMt
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt \
        -name Scer_S288c \
        --nuclear \
        > Scer_S288c.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Scer_S288c.seq.csv
    ```

## *Scer_wgs* WGS

1. Create `pop/scer_wgs.tsv` manually.

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Saccharomyces&order=organism
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid4930[Organism:exp]

    ```bash
    export GENUS_ID=4930
    export GENUS=saccharomyces
    mkdir -p ~/data/alignment/Fungi/scer_wgs          # operation directory
    mkdir -p ~/data/alignment/Fungi/GENOMES/scer_wgs  # sequence directory

    cd ~/data/alignment/Fungi/GENOMES/scer_wgs

    ... # paste codes from README.md

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

    Remove other species of Saccharomyces. Add three outgroup genomes.

    ```bash
    echo -e '#name\tprefix\torganism\tcontigs' > scer_wgs.tsv
    cat saccharomyces.tsv | grep Scer_ | sed "s/^Scer_//" >> scer_wgs.tsv
    rm saccharomyces.tsv

    # outgroups
    echo -e "Spar\tAABY\tSaccharomyces paradoxus NRRL Y-17217\t832" >> scer_wgs.tsv
    echo -e "Spas\tJTFI\tSaccharomyces pastorianus\t878" >> scer_wgs.tsv
    echo -e "Sbou\tJRHY\tSaccharomyces sp. boulardii\t193" >> scer_wgs.tsv
    ```

    Edit them to fix names and comment out bad strains.

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/scer_wgs
    cd ~/data/alignment/Fungi/GENOMES/scer_wgs

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/scer_wgs.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/scer_wgs.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download strains of *Saccharomyces cerevisiae* at good assembly status.

    Click the `Download table` link on the top-right of [Genome list](http://www.ncbi.nlm.nih.gov/genome/genomes/15),
    save it as .csv file.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/scer_wgs/DOWNLOAD
    cd ~/data/alignment/Fungi/GENOMES/scer_wgs/DOWNLOAD

    # Download S288c and EC1118 separately
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt \
        --nuclear -name S288c \
        > S288c.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000218975.1.assembly.txt \
        --nuclear --genbank --scaffold -name EC1118 \
        > EC1118.seq.csv

    echo "#strain_name,accession,strain_taxon_id,seq_name" > scer_wgs.seq.csv
    cat S288c.seq.csv EC1118.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> scer_wgs.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f scer_wgs.seq.csv

    ```

## *Scer_100* ASSEMBLY

1. Create working directory and download WGS strains as outgroups.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/scer_100
    cd ~/data/alignment/Fungi/GENOMES/scer_100

    echo -e '#name\tprefix\torganism\tcontigs' > scer_100.tsv

    # outgroups
    echo -e "Spar\tAABY\tSaccharomyces paradoxus NRRL Y-17217\t832" >> scer_100.tsv
    echo -e "Spas\tJTFI\tSaccharomyces pastorianus\t878" >> scer_100.tsv
    echo -e "Sbou\tJRHY\tSaccharomyces sp. boulardii\t193" >> scer_100.tsv

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f scer_100.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/scer_100.url.txt

    find WGS -name "*.gz" | xargs gzip -t

    ```

2. Download strains of *Saccharomyces cerevisiae* at good assembly status.

    Click the `Download table` link on the top-right of [Genome list](http://www.ncbi.nlm.nih.gov/genome/genomes/15),
    save it as .csv file.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/scer_100/DOWNLOAD
    cd ~/data/alignment/Fungi/GENOMES/scer_100/DOWNLOAD

    # Download S288c separately
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt \
        --nuclear -name S288c \
        > S288c.seq.csv

    mysql -ualignDB -palignDB ar_genbank -e "SELECT organism_name, species, assembly_accession FROM ar WHERE wgs_master = '' AND organism_name != species AND species_id = 4932" \
        | perl -nl -a -F"\t" -e '$n = $F[0]; $rx = quotemeta $F[1]; $n =~ s/$rx\s+//; $n =~ s/\W+/_/g; printf qq{%s\t%s\n}, $n, $F[2];' \
        | grep -v organism_name | grep -v S288c | grep -v EC1118 \
        | perl -nl -a -F"\t" -e '$str = q{echo } . $F[0] . qq{ \n}; $str .= q{perl ~/Scripts/withncbi/taxon/assembly_csv.pl} . qq{ \\\n}; $str .= q{-f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/} . $F[1] . qq{.assembly.txt \\\n}; $str .= q{--nuclear --genbank --chromosome -name } . $F[0] . qq{ \\\n}; $str .= q{>> non_wgs.seq.csv}; print $str . qq{\n}' \
        > ass_csv.sh

    echo > non_wgs.seq.csv
    sh ass_csv.sh

    echo "#strain_name,accession,strain_taxon_id,seq_name" > scer_100.seq.csv
    cat S288c.seq.csv non_wgs.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> scer_100.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f scer_100.seq.csv

    ```

## *Candida* WGS

1. Create `pop/candida.tsv` manually.

    Check NCBI pages

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=Candida*&order=organism
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid1535326[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid1535326[Organism:exp]

    The wgs page contains a lot of strains from other genus.

    ```bash
    export GENUS_ID=1535326
    export GENUS=candida
    mkdir -p ~/data/alignment/Fungi/$GENUS          # operation directory
    mkdir -p ~/data/alignment/Fungi/GENOMES/$GENUS  # sequence directory

    cd ~/data/alignment/Fungi/GENOMES/$GENUS

    ...
    # Edit raw2.tsv, remove lines containing CANDIDATUS or CANDIDATE DIVISION
    cat raw2.tsv | grep -v 'CANDIDATUS' | grep -v 'CANDIDATE DIVISION' > tmp.tsv
    mv tmp.tsv raw2.tsv
    ...

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/candida
    cd ~/data/alignment/Fungi/GENOMES/candida

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/candida.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/candida.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Pick targets.

    ```sql
    SELECT *
    FROM gr.gr
    WHERE genus_id = 1535326

    SELECT
        organism_name, assembly_accession
    FROM
        ar_refseq.ar
    WHERE
        genus_id = 1535326
    ORDER BY organism_name
    ```

    | assigned name  | organism_name                     | assembly_accession |
    | :------------  | :------------                     | :------------      |
    | Cdub_CD36      | *Candida dubliniensis* CD36       | GCF_000026945.1    |
    | Corh_Co_90_125 | *Candida orthopsilosis* Co 90-125 | GCF_000315875.1    |

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/candida/DOWNLOAD
    cd ~/data/alignment/Fungi/GENOMES/candida/DOWNLOAD

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000026945.1.assembly.txt \
        -name Cdub_CD36 \
        > Cdub_CD36.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000315875.1.assembly.txt \
        -name Corh_Co_90_125 \
        > Corh_Co_90_125.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Cdub_CD36.seq.csv

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Corh_Co_90_125.seq.csv
    ```

## *Fusarium* WGS

1. Create `pop/fusarium.tsv` manually.

    Check NCBI pages

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=fusarium*&order=organism
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid5506[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid5506[Organism:exp]

    This genus contains about 40 strains, use command lines listed in README.md.

    ```bash
    export GENUS_ID=5506
    export GENUS=fusarium
    mkdir -p ~/data/alignment/Fungi/$GENUS          # operation directory
    mkdir -p ~/data/alignment/Fungi/GENOMES/$GENUS  # sequence directory

    cd ~/data/alignment/Fungi/GENOMES/$GENUS

    ...

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/fusarium
    cd ~/data/alignment/Fungi/GENOMES/fusarium

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/fusarium.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/fusarium.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Pick targets.

    ```sql
    SELECT *
    FROM gr.gr
    WHERE genus_id = 5506
    ```

    | assigned name | organism_name                                | assembly_accession |
    | :------------ | :------------                                | :------------      |
    | Fgra_PH_1     | *Fusarium graminearum* PH-1                  | GCA_000240135.3    |
    | Foxy_4287     | *Fusarium oxysporum* f. sp. lycopersici 4287 | GCA_000149955.1    |
    | Fpse_CS3270   | *Fusarium pseudograminearum* CS3270          | GCA_000974265.1    |
    | Fver_7600     | *Fusarium verticillioides* 7600              | GCA_000149555.1    |

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/fusarium/DOWNLOAD
    cd ~/data/alignment/Fungi/GENOMES/fusarium/DOWNLOAD

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000240135.3.assembly.txt \
        -name Fgra_PH_1 \
        --genbank \
        --nuclear \
        > Fgra_PH_1.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000149955.1.assembly.txt \
        -name Foxy_4287 \
        --genbank \
        --nuclear \
        > Foxy_4287.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000974265.1.assembly.txt \
        -name Fpse_CS3270 \
        --genbank \
        --nuclear \
        > Fpse_CS3270.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000149555.1.assembly.txt \
        -name Fver_7600 \
        --genbank \
        --nuclear \
        > Fver_7600.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Fgra_PH_1.seq.csv

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Foxy_4287.seq.csv

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Fpse_CS3270.seq.csv

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Fver_7600.seq.csv
    ```

## *Aspergillus* WGS

1. Create `pop/aspergillus.tsv` manually.

    Check NCBI pages

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=aspergillus&order=organism
    * http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5052
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid5052[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid5052[Organism:exp]

    This genus contains about 40 strains, use command lines listed in README.md.

    ```bash
    export GENUS_ID=5052
    export GENUS=aspergillus
    mkdir -p ~/data/alignment/Fungi/$GENUS          # operation directory
    mkdir -p ~/data/alignment/Fungi/GENOMES/$GENUS  # sequence directory

    cd ~/data/alignment/Fungi/GENOMES/$GENUS

    ...

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/aspergillus
    cd ~/data/alignment/Fungi/GENOMES/aspergillus

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/aspergillus.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/aspergillus.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Pick targets.

    ```sql
    SELECT *
    FROM gr_euk.gr
    WHERE genus_id = 5052

    SELECT
        *
    FROM
        ar_genbank.ar
    WHERE
        genus_id = 5052
    ORDER BY organism_name
    ```

    | assigned name | organism_name                  | assembly_accession |
    | :------------ | :------------                  | :------------      |
    | Afum_Af293    | *Aspergillus fumigatus* Af293  | GCA_000002655.1    |
    | Anid_FGSC_A4  | *Aspergillus nidulans* FGSC A4 | GCA_000011425.1    |

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/aspergillus/DOWNLOAD
    cd ~/data/alignment/Fungi/GENOMES/aspergillus/DOWNLOAD

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000002655.1.assembly.txt \
        -name Afum_Af293 \
        --nuclear \
        > Afum_Af293.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000011425.1.assembly.txt \
        -name Anid_FGSC_A4 \
        --genbank \
        --nuclear \
        > Anid_FGSC_A4.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Afum_Af293.seq.csv

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Anid_FGSC_A4.seq.csv
    ```

## *Penicillium* WGS

1. Create `pop/penicillium.tsv` manually.

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=penicillium&order=organism
    * http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5073
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid5073[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid5073[Organism:exp]

    ```bash
    export GENUS_ID=5073
    export GENUS=penicillium
    mkdir -p ~/data/alignment/Fungi/$GENUS          # operation directory
    mkdir -p ~/data/alignment/Fungi/GENOMES/$GENUS  # sequence directory

    cd ~/data/alignment/Fungi/GENOMES/$GENUS

    ...

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Fungi/GENOMES/penicillium
    cd ~/data/alignment/Fungi/GENOMES/penicillium

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/penicillium.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/penicillium.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Pick targets.

    ```sql
    SELECT *
    FROM gr_euk.gr
    WHERE genus_id = 5073

    SELECT
        *
    FROM
        ar_genbank.ar
    WHERE
        genus_id = 5073
    ORDER BY organism_name
    ```

    There're no good target. Pchr_P2niaD18 is the only one on chromosome level, but is not de novo
    assembled and hasn't annotations.

    | assigned name | organism_name                      | assembly_accession |
    | :------------ | :------------                      | :------------      |
    | Pchr_P2niaD18 | *Penicillium chrysogenum* P2niaD18 | GCA_000710275.1    |

## *Plasmodium* WGS

    | name                      | taxon |
    | :---                      | :---  |
    | Plasmodium                | 5820  |
    | Plasmodium falciparum     | 5833  |
    | Plasmodium falciparum 3D7 | 36329 |

1. Create `pop/plasmodium.tsv` manually.

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=plasmodium&order=organism
    * http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5820
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid5820[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid5820[Organism:exp]

    ```bash
    export GENUS_ID=5820
    export GENUS=plasmodium
    mkdir -p ~/data/alignment/Protists/$GENUS          # operation directory
    mkdir -p ~/data/alignment/Protists/GENOMES/$GENUS  # sequence directory

    cd ~/data/alignment/Protists/GENOMES/$GENUS

    ...

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

    Remove all strains of Plasmodium falciparum. 3D7 will be injected later.

    ```bash
    mv plasmodium.tsv all.tsv
    cat all.tsv | grep -v Pfal_ > plasmodium.tsv
    ```

    Edit the tsv file to fix names and comment out bad strains.


2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Protists/GENOMES/plasmodium
    cd ~/data/alignment/Protists/GENOMES/plasmodium

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/plasmodium.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/plasmodium.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download *Plasmodium falciparum* 3D7.
  This step is totally manual operation. **Be careful.**

    | assigned name | organism_name               | assembly_accession           |
    | :------------ | :------------               | :------------                |
    | Pfal_3D7      | *Plasmodium falciparum* 3D7 | GCF_000002765.3.assembly.txt |

    ```bash
    mkdir -p ~/data/alignment/Protists/GENOMES/plasmodium/DOWNLOAD
    cd ~/data/alignment/Protists/GENOMES/plasmodium/DOWNLOAD

    # Omit chrMt
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000002765.3.assembly.txt \
        -name Pfal_3D7 \
        --nuclear \
        > Pfal_3D7.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Pfal_3D7.seq.csv
    ```

## *Plasmodium falciparum* WGS

1. Create `pop/Pfal.tsv` manually.

    ```bash
    export GENUS_ID=5820
    export GENUS=plasmodium
    mkdir -p ~/data/alignment/Protists/pfal          # operation directory
    mkdir -p ~/data/alignment/Protists/GENOMES/pfal  # sequence directory

    cd ~/data/alignment/Protists/GENOMES/pfal

    ... # paste codes from README.md

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

    Remove other species of Plasmodium. Add Prei as outgroup.

    ```bash
    echo -e '#name\tprefix\torganism\tcontigs' > pfal.tsv
    cat plasmodium.tsv | grep Pfal_ | sed "s/^Pfal_//" >> pfal.tsv
    rm plasmodium.tsv

    # outgroup
    echo -e "Prei\tCBXM\tPlasmodium reichenowi\t2055" >> pfal.tsv
    ```

    Edit the tsv file to fix names and comment out bad strains.

2. Create working directory and download WGS sequences.

    ```bash
    cd ~/data/alignment/Protists/GENOMES/pfal

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/pfal.tsv \
        --fix \
        --nofix Prei \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/pfal.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download *Plasmodium falciparum* 3D7.

    Click the `Download table` link on the top-right of [Genome list](http://www.ncbi.nlm.nih.gov/genome/genomes/33),
    save it as .csv file.

    ```bash
    mkdir -p ~/data/alignment/Protists/GENOMES/pfal/DOWNLOAD
    cd ~/data/alignment/Protists/GENOMES/pfal/DOWNLOAD

    # Download Plasmodium falciparum 3D7 separately (36329)
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000002765.3.assembly.txt \
        --nuclear -name 3D7 \
        > 3D7.seq.csv

    mysql -ualignDB -palignDB ar_genbank -e \
        "SELECT organism_name, species, assembly_accession FROM ar WHERE taxonomy_id IN (1036723, 5843, 1237627, 1036727, 1036725, 1036726, 57270, 5835, 1036724, 57266, 478859)" \
        | perl -nl -a -F"\t" -e '$n = $F[0]; $rx = quotemeta $F[1]; $n =~ s/$rx\s+//; $n =~ s/\W+/_/g; printf qq{%s\t%s\n}, $n, $F[2];' \
        | grep -v organism_name \
        | perl -nl -a -F"\t" -e '$str = q{echo } . $F[0] . qq{ \n}; $str .= q{perl ~/Scripts/withncbi/taxon/assembly_csv.pl} . qq{ \\\n}; $str .= q{-f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/} . $F[1] . qq{.assembly.txt \\\n}; $str .= q{ --scaffold --length 5000 --genbank -name } . $F[0] . qq{ \\\n}; $str .= q{>> non_wgs.seq.csv}; print $str . qq{\n}' \
        > ass_csv.sh

    echo > non_wgs.seq.csv
    sh ass_csv.sh

    echo "#strain_name,accession,strain_taxon_id,seq_name" > pfal.seq.csv
    cat 3D7.seq.csv non_wgs.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> pfal.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f pfal.seq.csv

    ```

## *Arabidopsis* 19 genomes

1. Sources.

    * [Project page](http://mus.well.ox.ac.uk/19genomes/)
    * [Download page](http://mus.well.ox.ac.uk/19genomes/fasta/)
    * [Paper](http://www.nature.com/nature/journal/v477/n7365/full/nature10414.html)

2. Download with my web page crawler.

    ```bash
    mkdir -p ~/data/alignment/others
    cd ~/data/alignment/others

    perl ~/Scripts/download/list.pl -u http://mus.well.ox.ac.uk/19genomes/fasta/
    perl ~/Scripts/download/download.pl -i 19genomes_fasta.yml -a
    aria2c -x 12 -s 4 -i /home/wangq/data/alignment/others/19genomes_fasta.yml.txt
    ```

3. *A. thaliana* and *A. lyrata* from ensembl genomes.

    Check [ensembl.md](ensembl.md). Ensembl data stored in `~/data/ensembl82`.

    ```bash
    # Atha
    mkdir -p ~/data/alignment/Ensembl/Atha
    cd ~/data/alignment/Ensembl/Atha

    find ~/data/ensembl82/fasta/arabidopsis_thaliana/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    mv Mt.fa Mt.fa.skip
    mv Pt.fa Pt.fa.skip

    # Alyr
    mkdir -p ~/data/alignment/Ensembl/Alyr
    cd ~/data/alignment/Ensembl/Alyr

    find ~/data/ensembl82/fasta/arabidopsis_lyrata/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-about toplevel.filtered.fa 10000000 .
    rm toplevel.fa toplevel.filtered.fa listFile
    ```

## *Orazy sativa* Japonica 24 genomes

1. Sources.

    * [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP003189)
    * [Paper](http://www.nature.com/nbt/journal/v30/n1/full/nbt.2050.html)

    Mapping strategy in [here](https://github.com/wang-q/sra/blob/master/japonica24_seq.pl).

2. 23 Japonica rices restore from previous .2bit files.

    I've used bwa-gatk pipeline to generate 23 Japonica rice genomes.

    Reference assembly of nipponbare in that time was MSU6, now is IRGSP-1.0. But I don't want to waste time to rebuild and RepeatMasker all sequences.

    ```bash
    find ~/data/alignment/rice/ -name "*.2bit" \
        | grep -v "_65" \
        | parallel basename {//} \
        | sort

    mkdir -p ~/data/alignment/others/japonica24
    cd ~/data/alignment/others/japonica24

    for d in IRGC11010 IRGC1107 IRGC12793 IRGC17757 IRGC2540 IRGC26872 IRGC27630 IRGC31856 IRGC32399 IRGC328 IRGC38698 IRGC38994 IRGC418 IRGC43325 IRGC43675 IRGC50448 IRGC55471 IRGC66756 IRGC8191 IRGC8244 IRGC9060 IRGC9062 RA4952
    do
        twoBitToFa ~/data/alignment/rice/${d}/chr.2bit ${d}.fa;
    done
    ```

3. nipponbare and 9311 from ensembl genomes.

    ```bash
    # OsatJap
    mkdir -p ~/data/alignment/Ensembl/OsatJap
    cd ~/data/alignment/Ensembl/OsatJap

    find ~/data/ensembl82/fasta/oryza_sativa/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    rm AC*.fa AP*.fa Syng*.fa

    # OsatInd
    mkdir -p ~/data/alignment/Ensembl/OsatInd
    cd ~/data/alignment/Ensembl/OsatInd

    find ~/data/ensembl82/fasta/oryza_indica/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    rm AA*.fa CH*.fa Sup*.fa
    ```

## *Drosophila* Population Genomics Project (dpgp)

1. Sources.

    * [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP005599)
    * [Paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003080)

    Mapping strategy in [here](https://github.com/wang-q/sra/blob/master/dpgp_seq.pl).

2. 21 genomes restore from previous .2bit files.

    ```bash
    find ~/data/alignment/dpgp/ -name "*.2bit" \
        | grep -v "_65" \
        | parallel basename {//} \
        | sort

    mkdir -p ~/data/alignment/others/dpgp
    cd ~/data/alignment/others/dpgp

    for d in CK1 CO15N CO2 CO4N ED10N EZ5N FR207 FR217 FR229 FR361 GA125 GA129 GA130 GA132 GA185 GU10 KN6 KR39 KT1 NG3N RC1 RG15 SP254 TZ8 UG7 UM526 ZI268 ZL130 ZO12 ZS37
    do
        twoBitToFa ~/data/alignment/dpgp/${d}/chr.2bit ${d}.fa;
    done
    ```

3. Dmel and Dsim from ensembl genomes.

    ```bash
    # Dmel
    mkdir -p ~/data/alignment/Ensembl/Dmel
    cd ~/data/alignment/Ensembl/Dmel

    find ~/data/ensembl82/fasta/drosophila_melanogaster/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    rm *Scaffold*.fa 211*.fa
    mv 4.fa 4.fa.skip
    mv Y.fa Y.fa.skip
    mv rDNA.fa rDNA.fa.skip
    mv dmel_mitochondrion_genome.fa dmel_mitochondrion_genome.fa.skip

    # Dsim
    mkdir -p ~/data/alignment/Ensembl/Dsim
    cd ~/data/alignment/Ensembl/Dsim

    find ~/data/ensembl82/fasta/drosophila_simulans/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-about toplevel.filtered.fa 10000000 .
    rm toplevel.fa toplevel.filtered.fa listFile
    ```

## Primates

1. Guild tree

    ```bash
    cd ~/data/alignment
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/hg19.100way.commonNames.nh
    tree_doctor hg19.100way.commonNames.nh --newick  \
        --prune-all-but \
        Human,Chimp,Gorilla,Orangutan,Gibbon,Rhesus,Crab_eating_macaque,Baboon,Green_monkey,Marmoset,Squirrel_monkey,Bushbaby,Chinese_tree_shrew \
        > primates_13way.nwk

    ```

2. All from ensembl.

    ```bash
    # Human
    mkdir -p ~/data/alignment/Ensembl/Human
    cd ~/data/alignment/Ensembl/Human

    find ~/data/ensembl82/fasta/homo_sapiens/dna/ -name "*dna_sm.primary_assembly*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    rm GL*.fa

    mv Y.fa Y.fa.skip
    mv MT.fa MT.fa.skip

    # Chimp
    mkdir -p ~/data/alignment/Ensembl/Chimp
    cd ~/data/alignment/Ensembl/Chimp

    find ~/data/ensembl82/fasta/pan_troglodytes/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat GL*.fa > Un.fa
    cat AACZ*.fa >> Un.fa

    rm GL*.fa AACZ*.fa

    # Gorilla
    mkdir -p ~/data/alignment/Ensembl/Gorilla
    cd ~/data/alignment/Ensembl/Gorilla

    find ~/data/ensembl82/fasta/gorilla_gorilla/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    # Orangutan
    mkdir -p ~/data/alignment/Ensembl/Orangutan
    cd ~/data/alignment/Ensembl/Orangutan

    find ~/data/ensembl82/fasta/pongo_abelii/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    # Rhesus
    mkdir -p ~/data/alignment/Ensembl/Rhesus
    cd ~/data/alignment/Ensembl/Rhesus

    find ~/data/ensembl82/fasta/macaca_mulatta/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat 1099*.fa > Un.fa

    rm 1099*.fa
    ```

## *Caenorhabditis elegans*

There are no suitable outgroups for *C. elegans*.

http://hgdownload.soe.ucsc.edu/goldenPath/ce10/multiz7way/ce10.commonNames.7way.nh

1. 40 wild strains from cele_mmp.

    Mapping strategy in [here](https://github.com/wang-q/sra/blob/master/cele_mmp_seq.pl).

	```bash
    mkdir -p ~/data/alignment/others/cele
    cd ~/data/alignment/others/cele

    find ~/data/dna-seq/cele_mmp/ -name "*.vcf.fasta" \
        | parallel -j 1 cp {} .
	```

2. Reference strain N2 from ensembl genomes

    ```bash
    # Dmel
    mkdir -p ~/data/alignment/Ensembl/Cele
    cd ~/data/alignment/Ensembl/Cele

    find ~/data/ensembl82/fasta/caenorhabditis_elegans/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
    faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    mv MtDNA.fa MtDNA.fa.skip
	```

## *Dictyostelium* WGS

| name                         | taxon  |
| :---                         | :---   |
| Dictyostelium                | 5782   |
| Dictyostelium discoideum     | 44689  |
| Dictyostelium discoideum AX4 | 352472 |

1. Create `pop/dictyostelium.tsv` manually.

    * http://www.ncbi.nlm.nih.gov/Traces/wgs/?page=1&term=dictyostelium&order=organism
    * http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=5782
    * http://www.ncbi.nlm.nih.gov/assembly?term=txid5782[Organism:exp]
    * http://www.ncbi.nlm.nih.gov/genome/?term=txid5782[Organism:exp]

    ```bash
    export GENUS_ID=5782
    export GENUS=dictyostelium
    mkdir -p ~/data/alignment/Protists/$GENUS          # operation directory
    mkdir -p ~/data/alignment/Protists/GENOMES/$GENUS  # sequence directory

    cd ~/data/alignment/Protists/GENOMES/$GENUS

    ...

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

    Remove Ddis AX4. AX4 will be injected later.

    ```bash
    mv dictyostelium.tsv all.tsv
    cat all.tsv | grep -v Ddis_ > dictyostelium.tsv
    ```

    Edit the tsv file to fix names and comment out bad strains.


2. Create working directory and download WGS sequences.

    ```bash
    mkdir -p ~/data/alignment/Protists/GENOMES/dictyostelium
    cd ~/data/alignment/Protists/GENOMES/dictyostelium

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f ~/Scripts/withncbi/pop/dictyostelium.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/dictyostelium.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download *Dictyostelium discoideum* AX4.
    This step is totally manual operation. **Be careful.**

| assigned name | organism_name                  | assembly_accession           |
| :------------ | :------------                  | :------------                |
| Ddis_AX4      | *Dictyostelium discoideum* AX4 | GCF_000004695.1.assembly.txt |

    ```bash
    mkdir -p ~/data/alignment/Protists/GENOMES/dictyostelium/DOWNLOAD
    cd ~/data/alignment/Protists/GENOMES/dictyostelium/DOWNLOAD

    # Omit MT and Ddp5 plasmid
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000004695.1.assembly.txt \
        -name Ddis_AX4 \
        --nuclear \
        | grep -v Ddp5 \
        > Ddis_AX4.seq.csv

    mysql -ualignDB -palignDB ar_genbank -e \
        "SELECT organism_name, species, assembly_accession FROM ar WHERE taxonomy_id IN (5786, 261658, 361076, 79012, 361072)" \
        | perl -nl -a -F"\t" -e '$n = $F[0]; $rx = quotemeta $F[1]; $n =~ s/$rx\s+//; $n =~ s/\W+/_/g; printf qq{%s\t%s\n}, $n, $F[2];' \
        | grep -v organism_name \
        | perl -nl -a -F"\t" -e '$str = q{echo } . $F[0] . qq{ \n}; $str .= q{perl ~/Scripts/withncbi/taxon/assembly_csv.pl} . qq{ \\\n}; $str .= q{-f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/} . $F[1] . qq{.assembly.txt \\\n}; $str .= q{ --scaffold --length 5000 --genbank -name } . $F[0] . qq{ \\\n}; $str .= q{>> non_wgs.seq.csv}; print $str . qq{\n}' \
        > ass_csv.sh

    echo > non_wgs.seq.csv
    sh ass_csv.sh

    echo "#strain_name,accession,strain_taxon_id,seq_name" > dictyostelium.seq.csv
    cat Ddis_AX4.seq.csv non_wgs.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> dictyostelium.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f dictyostelium.seq.csv

    ```

## *Dictyostelium discoideum*

1. Sources.

    * [SRA](http://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP002085)
    * [Reference genome](https://www.hgsc.bcm.edu/microbiome/dictyostelium-discoideum-ax4)

    Mapping strategy in [here](https://github.com/wang-q/sra/blob/master/dicty_seq.pl).

2. 18 genomes restore from previous .2bit files.

    One of which is AX4, the reference genome resequenced.

    ```bash
    find ~/data/alignment/dicty/ -name "*.2bit" \
        | grep -v "_65" \
        | parallel basename {//} \
        | sort

    mkdir -p ~/data/alignment/others/dicty
    cd ~/data/alignment/others/dicty

    for d in 68 70 AX4 QS11 QS17 QS18 QS23 QS36 QS37 QS4 QS69 QS73 QS74 QS80 QS9 S224 WS14 WS15
    do
        twoBitToFa ~/data/alignment/dicty/${d}/chr.2bit ${d}.fa;
    done
    ```

3. Ddis, Dfir, Dcit from NCBI.

    ```bash
    mkdir -p ~/data/alignment/Protists/GENOMES/Ddis/DOWNLOAD
    cd ~/data/alignment/Protists/GENOMES/Ddis/DOWNLOAD

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000004695.1.assembly.txt \
        -name AX4 \
        --nuclear \
        | grep -v Ddp5 \
        > AX4.seq.csv

    echo > non_wgs.seq.csv
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000286055.1.assembly.txt \
         --scaffold --length 5000 --genbank -name Dcit \
        >> non_wgs.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000277485.1.assembly.txt \
         --scaffold --length 5000 --genbank -name Dfir \
        >> non_wgs.seq.csv

    echo "#strain_name,accession,strain_taxon_id,seq_name" > Ddis.seq.csv
    cat AX4.seq.csv non_wgs.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> Ddis.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f Ddis.seq.csv

    ```
