<!-- TOC depth:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Operating steps for each groups](#operating-steps-for-each-groups)
	- [Align](#align)
		- [*Saccharomyces* WGS](#saccharomyces-wgs)
		- [*Candida* WGS](#candida-wgs)
		- [*Fusarium* WGS](#fusarium-wgs)
		- [*Aspergillus* WGS](#aspergillus-wgs)
		- [*Penicillium* WGS](#penicillium-wgs)
<!-- /TOC -->

# Operating steps for each groups

Less detailed than Trichoderma in [README.md](README.md), but include examples
for genomes out of WGS, which usually in better assembling levels.

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
