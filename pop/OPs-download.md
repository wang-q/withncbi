<!-- TOC depth:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Operating steps for each groups](#operating-steps-for-each-groups)
	- [Download](#download)
		- [*Saccharomyces* WGS](#saccharomyces-wgs)
		- [*Scer_wgs* WGS](#scerwgs-wgs)
		- [*Scer_100* ASSEMBLY](#scer100-assembly)
		- [*Candida* WGS](#candida-wgs)
		- [*Fusarium* WGS](#fusarium-wgs)
		- [*Aspergillus* WGS](#aspergillus-wgs)
		- [*Penicillium* WGS](#penicillium-wgs)
<!-- /TOC -->

# Operating steps for each groups

Less detailed than Trichoderma in [README.md](README.md), but include examples
for genomes out of WGS, which usually in better assembling levels.

## Download

### *Saccharomyces* WGS

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

### *Scer_wgs* WGS

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

### *Scer_100* ASSEMBLY

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

### *Candida* WGS

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

### *Fusarium* WGS

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

### *Aspergillus* WGS

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

### *Penicillium* WGS

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
