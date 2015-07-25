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

    ...

    # Cleaning
    rm raw*.*sv
    unset GENUS_ID
    unset GENUS
    ```

    Remove all strains of Saccharomyces cerevisiae. S288c will be injected later.

    ```bash
    mv saccharomyces.tsv all.tsv
    echo -e '#name\tprefix\torganism\tcontigs' > scer_new.tsv
    cat all.tsv | grep Scer_ >> scer_new.tsv
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

    aria2c -x 6 -s 3 -c -i WGS/saccharomyces.url.txt

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
        -f Scer_S288c.seq.csv  \
        -r -p
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

    aria2c -x 6 -s 3 -c -i WGS/candida.url.txt

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
        -f Cdub_CD36.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -f Corh_Co_90_125.seq.csv  \
        -r -p
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

    aria2c -x 6 -s 3 -c -i WGS/fusarium.url.txt

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
        -f Fgra_PH_1.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -f Foxy_4287.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -f Fpse_CS3270.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -f Fver_7600.seq.csv  \
        -r -p
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

    aria2c -x 6 -s 3 -c -i WGS/aspergillus.url.txt

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
        -f Afum_Af293.seq.csv  \
        -r -p

    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -f Anid_FGSC_A4.seq.csv  \
        -r -p
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

    aria2c -x 6 -s 3 -c -i WGS/penicillium.url.txt

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

## Align

### *Saccharomyces* WGS

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/Fungi/saccharomyces
    cd ~/data/alignment/Fungi/saccharomyces

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/GENOMES/saccharomyces/WGS/saccharomyces.data.yml \
        -o ~/Scripts/withncbi/pop/saccharomyces_test.yml \
        -d ~/data/alignment/Fungi/GENOMES/saccharomyces/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=saccharomyces \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir='~/data/alignment/Fungi/saccharomyces' \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/Fungi/GENOMES/saccharomyces/DOWNLOAD \
        --downloaded 'name=Scer_S288c;taxon=559292;sciname=Saccharomyces cerevisiae S288c' \
        --plan 'name=four_way;t=Scer_S288c;qs=Sbou_ATCC_MYA_796,Spar_NRRL_Y_17217,Spas_CBS_1483'
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/saccharomyces_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_four_way.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

### *Candida* WGS

1. `gen_pop_conf.pl`

    Pay attentions to --downloaded orders. The first one will be the default target.

    ```bash
    mkdir -p ~/data/alignment/Fungi/candida
    cd ~/data/alignment/Fungi/candida

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/candida/WGS/candida.data.yml \
        -o ~/Scripts/withncbi/pop/candida_test.yml \
        -d ~/data/alignment/Fungi/candida/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=candida \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir='~/data/alignment/Fungi/candida' \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/Fungi/candida/DOWNLOAD \
        --downloaded 'name=Cdub_CD36;taxon=573826;sciname=Candida dubliniensis CD36' \
        --downloaded 'name=Corh_Co_90_125;taxon=1136231;sciname=Candida orthopsilosis Co 90-125' \
        --plan 'name=four_way;t=Cdub_CD36;qs=Corh_Co_90_125,Calb_WO_1,Ctro_MYA_3404' \
        --plan 'name=four_way_2;t=Corh_Co_90_125;qs=Cdub_CD36,Calb_WO_1,Ctro_MYA_3404'
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/candida_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_four_way.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_four_way_2.sh

    sh 3_pair_cmd.sh # Only do this when target switched
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

### *Fusarium* WGS

1. `gen_pop_conf.pl`

    Pay attentions to --downloaded orders. The first one will be the default target.

    ```bash
    mkdir -p ~/data/alignment/Fungi/fusarium
    cd ~/data/alignment/Fungi/fusarium

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/fusarium/WGS/fusarium.data.yml \
        -o ~/Scripts/withncbi/pop/fusarium_test.yml \
        -d ~/data/alignment/Fungi/fusarium/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=fusarium \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir='~/data/alignment/Fungi/fusarium' \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/Fungi/fusarium/DOWNLOAD \
        --downloaded 'name=Cdub_CD36;taxon=573826;sciname=Candida dubliniensis CD36' \
        --downloaded 'name=Corh_Co_90_125;taxon=1136231;sciname=Candida orthopsilosis Co 90-125' \
        --plan 'name=four_way;t=Cdub_CD36;qs=Corh_Co_90_125,Calb_WO_1,Ctro_MYA_3404' \
        --plan 'name=four_way_2;t=Corh_Co_90_125;qs=Cdub_CD36,Calb_WO_1,Ctro_MYA_3404'
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/fusarium_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_XXX.sh

    # sh 3_pair_cmd.sh # Only do this when target switched, e.g. four_way_2
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

### *Aspergillus* WGS

1. `gen_pop_conf.pl`

    Pay attentions to --downloaded orders. The first one will be the default target.

    ```bash
    mkdir -p ~/data/alignment/Fungi/aspergillus
    cd ~/data/alignment/Fungi/aspergillus

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/aspergillus/WGS/aspergillus.data.yml \
        -o ~/Scripts/withncbi/pop/aspergillus_test.yml \
        -d ~/data/alignment/Fungi/aspergillus/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=aspergillus \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir='~/data/alignment/Fungi/aspergillus' \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/Fungi/aspergillus/DOWNLOAD \
        --downloaded 'name=Afum_Af293;taxon=330879;sciname=Aspergillus fumigatus Af293' \
        --downloaded 'name=Anid_FGSC_A4;taxon=227321;sciname=Aspergillus nidulans FGSC A4' \
        --plan 'name=Afum_7way;t=Afum_Af293;qs=Afum_A1163,Afum_AF10,Afum_AF210,Afum_Af293,Afum_niveus,Afum_Z5' \
        --plan 'name=four_way_2;t=Anid_FGSC_A4;qs=Afum_Af293,Anig_ATCC_1015,Aory_RIB40'
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/aspergillus_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Afum_7way.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_four_way_2.sh

    sh 3_pair_cmd.sh # Only do this when target switched
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

### *Penicillium* WGS

1. `gen_pop_conf.pl`

    Pay attentions to --downloaded orders. The first one will be the default target.

    ```bash
    mkdir -p ~/data/alignment/Fungi/penicillium
    cd ~/data/alignment/Fungi/penicillium

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/penicillium/WGS/penicillium.data.yml \
        -o ~/Scripts/withncbi/pop/penicillium_test.yml \
        -d ~/data/alignment/Fungi/penicillium/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=penicillium \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir='~/data/alignment/Fungi/penicillium' \
        --opt rm_species=Fungi \
        --opt per_seq_min_contig=30000 \
        --per_seq Pchr_P2niaD18
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/penicillium_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```
