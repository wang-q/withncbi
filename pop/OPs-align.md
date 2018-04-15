# Aligning steps for each groups

Less detailed than Trichoderma in [README.md](README.md), but include examples for genomes out of
WGS, which usually in better assembling levels.

[TOC levels=1-3]: # " "
- [Aligning steps for each groups](#aligning-steps-for-each-groups)
- [*Saccharomyces* WGS](#saccharomyces-wgs)
- [Scer_wgs WGS](#scer_wgs-wgs)
- [Scer_100 ASSEMBLY](#scer_100-assembly)
- [*Candida* WGS](#candida-wgs)
- [*Fusarium* WGS](#fusarium-wgs)
- [*Aspergillus* WGS](#aspergillus-wgs)
- [*Penicillium* WGS](#penicillium-wgs)
- [*Plasmodium* WGS](#plasmodium-wgs)
- [*Plasmodium falciparum* WGS](#plasmodium-falciparum-wgs)
- [*Arabidopsis* 19 genomes](#arabidopsis-19-genomes)
- [*Orazy sativa* Japonica 24 genomes](#orazy-sativa-japonica-24-genomes)
- [*Drosophila* Population Genomics Project (dpgp)](#drosophila-population-genomics-project-dpgp)
- [Primates](#primates)
- [Human individuals from Simons project](#human-individuals-from-simons-project)
- [*Caenorhabditis elegans* million mutation project (cele_mmp)](#caenorhabditis-elegans-million-mutation-project-cele_mmp)
- [*Dictyostelium* WGS](#dictyostelium-wgs)
- [*Dictyostelium discoideum*](#dictyostelium-discoideum)
- [Mouse](#mouse)


# *Saccharomyces* WGS

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
        --download 'name=Scer_S288c;taxon=559292;sciname=Saccharomyces cerevisiae S288c' \
        --download 'name=Seub_FM1318;taxon=1080349;sciname=Saccharomyces eubayanus FM1318' \
        --plan 'name=plan_test;t=Scer_S288c;qs=Spar_NRRL_Y_17217,Spas_CBS_1483,Ssp_ATCC_MYA_796,Seub_FM1318'

    ```

2. Rest routing things.

    ```bash
    cd ~/data/alignment/Fungi/saccharomyces

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i ~/Scripts/withncbi/pop/saccharomyces_test.yml

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
    sh plan_plan_test.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

3. Create a summary xlsx.

Manually combine `~/data/alignment/Fungi/GENOMES/saccharomyces/WGS/saccharomyces.csv` and
`~/data/alignment/Fungi/saccharomyces/basicstat.xlsx`.

# Scer_wgs WGS

1. `gen_pop_conf.pl`

    ```bash
    # create downloaded genome list
    cat ~/data/alignment/Fungi/GENOMES/scer_wgs/DOWNLOAD/scer_wgs.seq.csv \
        | grep -v "^#" \
        | cut -d',' -f1,3 \
        | uniq \
        | perl -nl -a -F"," -e 'printf qq{    --download "name=%s;taxon=%s" \\\n}, $F[0], $F[1];'

    mkdir -p ~/data/alignment/Fungi/scer_wgs
    cd ~/data/alignment/Fungi/scer_wgs

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/GENOMES/scer_wgs/WGS/scer_wgs.data.yml \
        -o ~/Scripts/withncbi/pop/scer_wgs_test.yml \
        -d ~/data/alignment/Fungi/GENOMES/scer_wgs/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=scer_wgs \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir="~/data/alignment/Fungi/scer_wgs" \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/Fungi/GENOMES/scer_wgs/DOWNLOAD \
        --download "name=S288c;taxon=559292" \
        --download "name=RM11_1a;taxon=285006" \
        --download "name=EC1118;taxon=643680" \
        --plan 'name=Scer_n7_pop;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789' \
        --plan 'name=Scer_n7_Spar;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789,Spar;o=Spar' \
        --plan 'name=Scer_n7_Sbou;t=S288c;qs=EC1118,Kyokai_no_7,RM11_1a,Sigma1278b,T7,YJM789,Sbou;o=Sbou' \
        -y
    ```

2. Rest routing things.

    ```bash
    cd ~/data/alignment/Fungi/scer_wgs

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/scer_wgs_test.yml

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
    sh plan_Scer_n7_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Scer_n7_Spar.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Scer_n7_Sbou.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

3. Create a summary xlsx.

Manually combine `~/data/alignment/Fungi/GENOMES/scer_wgs/WGS/scer_wgs.csv` and
`~/data/alignment/Fungi/scer_wgs/basicstat.xlsx`.

# Scer_100 ASSEMBLY

1. `gen_pop_conf.pl`

    ```bash
    export GROUP_NAME=scer_100

    # create downloaded genome list
    cat ~/data/alignment/Fungi/GENOMES/${GROUP_NAME}/DOWNLOAD/${GROUP_NAME}.seq.csv \
        | grep -v "^#" \
        | cut -d',' -f1,3 \
        | uniq \
        | perl -nl -a -F"," -e 'printf qq{    --download "name=%s;taxon=%s" \\\n}, $F[0], $F[1];'

    cat ~/data/alignment/Fungi/GENOMES/${GROUP_NAME}/DOWNLOAD/${GROUP_NAME}.seq.csv \
        | grep -v "^#" \
        | cut -d',' -f1,3 \
        | uniq \
        | perl -nl -a -F"," -e 'printf qq{%s,}, $F[0];'

    mkdir -p ~/data/alignment/Fungi/${GROUP_NAME}
    cd ~/data/alignment/Fungi/${GROUP_NAME}

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/GENOMES/${GROUP_NAME}/WGS/${GROUP_NAME}.data.yml \
        -o ~/Scripts/withncbi/pop/${GROUP_NAME}_test.yml \
        -d ~/data/alignment/Fungi/GENOMES/${GROUP_NAME}/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=${GROUP_NAME} \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir="~/data/alignment/Fungi/${GROUP_NAME}" \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/Fungi/GENOMES/${GROUP_NAME}/DOWNLOAD \
        --download "name=S288c;taxon=559292" \
        --download "name=YJM993;taxon=1294331" \
        --download "name=YJM1078;taxon=1296266" \
        --download "name=YJM195;taxon=1294305" \
        --download "name=YJM270;taxon=1294308" \
        --download "name=YJM470;taxon=1294313" \
        --download "name=YJM683;taxon=1294320" \
        --download "name=YJM689;taxon=1294321" \
        --download "name=YJM693;taxon=1294322" \
        --download "name=YJM1248;taxon=1294340" \
        --download "name=YJM1252;taxon=1294342" \
        --download "name=YJM1273;taxon=1294343" \
        --download "name=YJM1342;taxon=1294352" \
        --download "name=YJM1385;taxon=1294357" \
        --download "name=YJM1387;taxon=1294359" \
        --download "name=YJM1388;taxon=1294360" \
        --download "name=YJM1389;taxon=1294361" \
        --download "name=YJM1399;taxon=1294362" \
        --download "name=YJM1402;taxon=1294365" \
        --download "name=YJM1418;taxon=1294368" \
        --download "name=YJM1439;taxon=1294372" \
        --download "name=YJM1443;taxon=1294373" \
        --download "name=YJM1444;taxon=1294374" \
        --download "name=YJM1447;taxon=1294375" \
        --download "name=YJM1460;taxon=1294377" \
        --download "name=YJM1549;taxon=1294384" \
        --download "name=YJM1573;taxon=1294385" \
        --download "name=YJM1592;taxon=1294387" \
        --download "name=YJM244;taxon=1294306" \
        --download "name=YJM1083;taxon=1292971" \
        --download "name=YJM1129;taxon=1293430" \
        --download "name=YJM189;taxon=1294303" \
        --download "name=YJM193;taxon=1294304" \
        --download "name=YJM248;taxon=1294307" \
        --download "name=YJM271;taxon=1294309" \
        --download "name=YJM320;taxon=947042" \
        --download "name=YJM326;taxon=468558" \
        --download "name=YJM428;taxon=947044" \
        --download "name=YJM450;taxon=1294310" \
        --download "name=YJM451;taxon=502869" \
        --download "name=YJM453;taxon=1294311" \
        --download "name=YJM456;taxon=1294312" \
        --download "name=YJM541;taxon=1294314" \
        --download "name=YJM554;taxon=1294315" \
        --download "name=YJM555;taxon=1294316" \
        --download "name=YJM627;taxon=1294317" \
        --download "name=YJM681;taxon=1294318" \
        --download "name=YJM682;taxon=1294319" \
        --download "name=YJM969;taxon=1294323" \
        --download "name=YJM972;taxon=1294324" \
        --download "name=YJM975;taxon=1294325" \
        --download "name=YJM978;taxon=1294326" \
        --download "name=YJM981;taxon=1294327" \
        --download "name=YJM984;taxon=1294328" \
        --download "name=YJM987;taxon=1294329" \
        --download "name=YJM990;taxon=1294330" \
        --download "name=YJM996;taxon=1294332" \
        --download "name=YJM1133;taxon=1294333" \
        --download "name=YJM1190;taxon=1294334" \
        --download "name=YJM1199;taxon=1294335" \
        --download "name=YJM1202;taxon=1294336" \
        --download "name=YJM1208;taxon=1294337" \
        --download "name=YJM1242;taxon=1294338" \
        --download "name=YJM1244;taxon=1294339" \
        --download "name=YJM1250;taxon=1294341" \
        --download "name=YJM1307;taxon=1294345" \
        --download "name=YJM1311;taxon=1294346" \
        --download "name=YJM1326;taxon=1294347" \
        --download "name=YJM1332;taxon=1294348" \
        --download "name=YJM1336;taxon=1294349" \
        --download "name=YJM1338;taxon=1294350" \
        --download "name=YJM1341;taxon=1294351" \
        --download "name=YJM1355;taxon=1294353" \
        --download "name=YJM1356;taxon=1294354" \
        --download "name=YJM1381;taxon=1294355" \
        --download "name=YJM1383;taxon=1294356" \
        --download "name=YJM1386;taxon=1294358" \
        --download "name=YJM1400;taxon=1294363" \
        --download "name=YJM1401;taxon=1294364" \
        --download "name=YJM1415;taxon=1294366" \
        --download "name=YJM1417;taxon=1294367" \
        --download "name=YJM1419;taxon=1294369" \
        --download "name=YJM1433;taxon=1294370" \
        --download "name=YJM1450;taxon=1294376" \
        --download "name=YJM1463;taxon=1294378" \
        --download "name=YJM1477;taxon=1294379" \
        --download "name=YJM1478;taxon=1294380" \
        --download "name=YJM1479;taxon=1294381" \
        --download "name=YJM1526;taxon=1294382" \
        --download "name=YJM1527;taxon=1294383" \
        --download "name=YJM1574;taxon=1294386" \
        --download "name=YJM1615;taxon=1294388" \
        --download "name=YJM1304;taxon=1294344" \
        --download "name=YJM1434;taxon=1294371" \
        --plan 'name=plan_og_spar;t=S288c;qs=YJM993,YJM1078,YJM195,YJM270,YJM470,YJM683,YJM689,YJM693,YJM1248,YJM1252,YJM1273,YJM1342,YJM1385,YJM1387,YJM1388,YJM1389,YJM1399,YJM1402,YJM1418,YJM1439,YJM1443,YJM1444,YJM1447,YJM1460,YJM1549,YJM1573,YJM1592,YJM244,YJM1083,YJM1129,YJM189,YJM193,YJM248,YJM271,YJM320,YJM326,YJM428,YJM450,YJM451,YJM453,YJM456,YJM541,YJM554,YJM555,YJM627,YJM681,YJM682,YJM969,YJM972,YJM975,YJM978,YJM981,YJM984,YJM987,YJM990,YJM996,YJM1133,YJM1190,YJM1199,YJM1202,YJM1208,YJM1242,YJM1244,YJM1250,YJM1307,YJM1311,YJM1326,YJM1332,YJM1336,YJM1338,YJM1341,YJM1355,YJM1356,YJM1381,YJM1383,YJM1386,YJM1400,YJM1401,YJM1415,YJM1417,YJM1419,YJM1433,YJM1450,YJM1463,YJM1477,YJM1478,YJM1479,YJM1526,YJM1527,YJM1574,YJM1615;o=Spar'

    unset GROUP_NAME
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i ~/Scripts/withncbi/pop/scer_100_test.yml

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
    sh plan_og_spar.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

# *Candida* WGS

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
        --download 'name=Cdub_CD36;taxon=573826;sciname=Candida dubliniensis CD36' \
        --download 'name=Corh_Co_90_125;taxon=1136231;sciname=Candida orthopsilosis Co 90-125' \
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

# *Fusarium* WGS

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
        --download 'name=Cdub_CD36;taxon=573826;sciname=Candida dubliniensis CD36' \
        --download 'name=Corh_Co_90_125;taxon=1136231;sciname=Candida orthopsilosis Co 90-125' \
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

# *Aspergillus* WGS

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
        --download 'name=Afum_Af293;taxon=330879;sciname=Aspergillus fumigatus Af293' \
        --download 'name=Anid_FGSC_A4;taxon=227321;sciname=Aspergillus nidulans FGSC A4' \
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

# *Penicillium* WGS

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

# *Plasmodium* WGS

1. RM species

    It't OK to not specifies RM species. Protists have very few repeats records.

    ```bash
    /usr/local/Cellar/repeatmasker/4.0.5/libexec/util/queryTaxonomyDatabase.pl -taxDBFile /usr/local/Cellar/repeatmasker/4.0.5/libexec/Libraries/taxonomy.dat -species Apicomplexa
    /usr/local/Cellar/repeatmasker/4.0.5/libexec/util/queryRepeatDatabase.pl -species Apicomplexa -stat
    ```

2. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/Protists/plasmodium
    cd ~/data/alignment/Protists/plasmodium

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Protists/GENOMES/plasmodium/WGS/plasmodium.data.yml \
        -o ~/Scripts/withncbi/pop/plasmodium_test.yml \
        -d ~/data/alignment/Protists/GENOMES/plasmodium/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=plasmodium \
        --opt base_dir='~/data/alignment/Protists' \
        --opt data_dir='~/data/alignment/Protists/plasmodium' \
        --dd ~/data/alignment/Protists/GENOMES/plasmodium/DOWNLOAD \
        --download 'name=Pfal_3D7;taxon=36329;sciname=Plasmodium falciparum 3D7' \
        -y
    ```

3. Rest routing things.

    ```bash
    cd ~/data/alignment/Protists/plasmodium

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i ~/Scripts/withncbi/pop/plasmodium_test.yml

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

4. Pick outgroup.

    Prei is the only one.

# *Plasmodium falciparum* WGS

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/Protists/pfal
    cd ~/data/alignment/Protists/pfal

    cat ~/data/alignment/Protists/GENOMES/pfal/DOWNLOAD/pfal.seq.csv \
        | grep -v "^#" \
        | cut -d',' -f1,3 \
        | uniq \
        | perl -nl -a -F"," -e 'printf qq{    --download "name=%s;taxon=%s" \\\n}, $F[0], $F[1];'

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Protists/GENOMES/pfal/WGS/pfal.data.yml \
        -o ~/Scripts/withncbi/pop/pfal_test.yml \
        -d ~/data/alignment/Protists/GENOMES/pfal/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=pfal \
        --opt base_dir='~/data/alignment/Protists' \
        --opt data_dir='~/data/alignment/Protists/pfal' \
        --dd ~/data/alignment/Protists/GENOMES/pfal/DOWNLOAD \
        --download 'name=3D7;taxon=36329;sciname=Plasmodium falciparum 3D7' \
        --download "name=CAMP_Malaysia;taxon=5835" \
        --download "name=NF54;taxon=5843" \
        --download "name=7G8;taxon=57266" \
        --download "name=Palo_Alto_Uganda;taxon=57270" \
        --download "name=Santa_Lucia;taxon=478859" \
        --download "name=Vietnam_Oak_Knoll_FVO_;taxon=1036723" \
        --download "name=FCH_4;taxon=1036724" \
        --download "name=Tanzania_2000708_;taxon=1036725" \
        --download "name=NF135_5_C10;taxon=1036726" \
        --download "name=MaliPS096_E11;taxon=1036727" \
        --download "name=UGT5_1;taxon=1237627" \
        --plan 'name=Pfal_n15_Prei;t=3D7;qs=7G8,CAMP_Malaysia,Dd2,FCH_4,HB3,IGH_CR14,MaliPS096_E11,NF135_5_C10,NF54,Palo_Alto_Uganda,Santa_Lucia,Tanzania_2000708_,UGT5_1,Vietnam_Oak_Knoll_FVO_,Prei;o=Prei' \
        --plan 'name=Pfal_n11_Prei;t=3D7;qs=7G8,CAMP_Malaysia,HB3,IGH_CR14,MaliPS096_E11,NF135_5_C10,Palo_Alto_Uganda,Santa_Lucia,Tanzania_2000708_,UGT5_1,Vietnam_Oak_Knoll_FVO_,Prei;o=Prei' \
        --plan 'name=Pfal_n11_pop;t=3D7;qs=7G8,CAMP_Malaysia,HB3,IGH_CR14,MaliPS096_E11,NF135_5_C10,Palo_Alto_Uganda,Santa_Lucia,Tanzania_2000708_,UGT5_1,Vietnam_Oak_Knoll_FVO_' \
        --plan 'name=Pfal_n7_pop;t=3D7;qs=7G8,CAMP_Malaysia,HB3,NF135_5_C10,Santa_Lucia,Vietnam_Oak_Knoll_FVO_' \
        --plan 'name=Pfal_n7_Prei;t=3D7;qs=7G8,CAMP_Malaysia,HB3,NF135_5_C10,Santa_Lucia,Vietnam_Oak_Knoll_FVO_,Prei;o=Prei' \
        -y
    ```

2. Rest routing things.

    ```bash
    cd ~/data/alignment/Protists/pfal

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i ~/Scripts/withncbi/pop/pfal_test.yml

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
    sh plan_Pfal_n15_Prei.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Pfal_n11_Prei.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Pfal_n11_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Pfal_n7_Prei.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Pfal_n7_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

# *Arabidopsis* 19 genomes

1. Create data.yml manually.

```bash
mkdir -p ~/data/alignment/arabidopsis82
cd ~/data/alignment/arabidopsis82

cat <<EOF > arabidopsis82_data.yml
---
data:
  - coverage: 25
    name: Bur_0
    origin: Ireland
    original_id: 3702
  - coverage: 47
    name: Can_0
    origin: Canary Isles
    original_id: 3702
  - coverage: 50
    name: Ct_1
    origin: Italy
    original_id: 3702
  - coverage: 52
    name: Edi_0
    origin: Scotland
    original_id: 3702
  - coverage: 33
    name: Hi_0
    origin: Netherlands
    original_id: 3702
  - coverage: 28
    name: Kn_0
    origin: Lithuania
    original_id: 3702
  - coverage: 27
    name: Ler_0
    origin: Poland
    original_id: 3702
  - coverage: 30
    name: Mt_0
    origin: Libya
    original_id: 3702
  - coverage: 38
    name: No_0
    origin: Germany
    original_id: 3702
  - coverage: 54
    name: Oy_0
    origin: Norway
    original_id: 3702
  - coverage: 41
    name: Po_0
    origin: Germany
    original_id: 3702
  - coverage: 38
    name: Rsch_4
    origin: Russia
    original_id: 3702
  - coverage: 40
    name: Sf_2
    origin: Spain
    original_id: 3702
  - coverage: 48
    name: Tsu_0
    origin: Japan
    original_id: 3702
  - coverage: 40
    name: Wil_2
    origin: Russia
    original_id: 3702
  - coverage: 33
    name: Ws_0
    origin: Russia
    original_id: 3702
  - coverage: 26
    name: Wu_0
    origin: Germany
    original_id: 3702
  - coverage: 31
    name: Zu_0
    origin: Germany
    original_id: 3702
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/arabidopsis82
    cd ~/data/alignment/arabidopsis82

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i arabidopsis82_data.yml \
        -o ~/Scripts/withncbi/pop/arabidopsis82_test.yml \
        -d ~/data/alignment/others/19genomes/fasta/MASKED \
        -m name \
        -r '*.fas' \
        --opt group_name=arabidopsis82 \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/arabidopsis82' \
        --dd ~/data/alignment/Ensembl \
        --download 'name=Atha;taxon=3702;sciname=Arabidopsis thaliana' \
        --download 'name=Alyr;taxon=59689;sciname=Arabidopsis lyrata' \
        --plan 'name=Ath_n19_pop;t=Atha;qs=Bur_0,Can_0,Ct_1,Edi_0,Hi_0,Kn_0,Ler_0,Mt_0,No_0,Oy_0,Po_0,Rsch_4,Sf_2,Tsu_0,Wil_2,Ws_0,Wu_0,Zu_0' \
        --plan 'name=Ath_n19_Alyr;t=Atha;qs=Bur_0,Can_0,Ct_1,Edi_0,Hi_0,Kn_0,Ler_0,Mt_0,No_0,Oy_0,Po_0,Rsch_4,Sf_2,Tsu_0,Wil_2,Ws_0,Wu_0,Zu_0,Alyr;o=Alyr' \
        -y
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/arabidopsis82_test.yml

    sh 01_file.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Ath_n19_pop.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Ath_n19_Alyr.sh

    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

# *Orazy sativa* Japonica 24 genomes

1. Create data.yml manually.

```bash
mkdir -p ~/data/alignment/rice82
cd ~/data/alignment/rice82

cat <<EOF > rice82_data.yml
---
data:
  - coverage: 11.33
    name: IRGC1107
    original_id: 4530
    subgroup: TEJ
  - coverage: 17.49
    name: IRGC2540
    original_id: 4530
    subgroup: TEJ
  - coverage: 9.6
    name: IRGC27630
    original_id: 4530
    subgroup: TEJ
  - coverage: 13.21
    name: IRGC32399
    original_id: 4530
    subgroup: TEJ
  - coverage: 13.27
    name: IRGC418
    original_id: 4530
    subgroup: TEJ
  - coverage: 15.99
    name: IRGC55471
    original_id: 4530
    subgroup: TEJ
  - coverage: 13.75
    name: IRGC8191
    original_id: 4530
    subgroup: TEJ
  - coverage: 16.07
    name: IRGC38698
    original_id: 4530
    subgroup: TEJ
  - coverage: 11.51
    name: IRGC11010
    original_id: 4530
    subgroup: TRJ
  - coverage: 10.74
    name: IRGC17757
    original_id: 4530
    subgroup: TRJ
  - coverage: 12.5
    name: IRGC328
    original_id: 4530
    subgroup: TRJ
  - coverage: 11.79
    name: IRGC43325
    original_id: 4530
    subgroup: TRJ
  - coverage: 10.97
    name: IRGC43675
    original_id: 4530
    subgroup: TRJ
  - coverage: 15.35
    name: IRGC50448
    original_id: 4530
    subgroup: TRJ
  - coverage: 11.83
    name: IRGC66756
    original_id: 4530
    subgroup: TRJ
  - coverage: 11.27
    name: IRGC8244
    original_id: 4530
    subgroup: TRJ
  - coverage: 12.78
    name: IRGC26872
    original_id: 4530
    subgroup: TRJ
  - coverage: 14.64
    name: IRGC12793
    original_id: 4530
    subgroup: ARO
  - coverage: 13.6
    name: IRGC38994
    original_id: 4530
    subgroup: ARO
  - coverage: 12.13
    name: IRGC9060
    original_id: 4530
    subgroup: ARO
  - coverage: 13.11
    name: IRGC9062
    original_id: 4530
    subgroup: ARO
  - coverage: 13.98
    name: RA4952
    original_id: 4530
    subgroup: ARO
  - coverage: 12.45
    name: IRGC31856
    original_id: 4530
    subgroup: ARO
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/rice82
    cd ~/data/alignment/rice82

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i rice82_data.yml \
        -o ~/Scripts/withncbi/pop/rice82_test.yml \
        -d ~/data/alignment/others/japonica24 \
        -m name \
        -r '*.fa' \
        --opt group_name=rice82 \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/rice82' \
        --dd ~/data/alignment/Ensembl \
        --download 'name=OsatJap;taxon=39947;sciname=Oryza sativa Japonica' \
        --download 'name=OsatInd;taxon=39946;sciname=Oryza sativa Indica' \
        --plan 'name=OsatJap_n8_OsatInd;t=OsatJap;qs=IRGC2540,IRGC38698,IRGC50448,IRGC26872,IRGC12793,RA4952,OsatInd;o=OsatInd' \
        --plan 'name=OsatJap_n24_pop;t=OsatJap;qs=IRGC1107,IRGC2540,IRGC27630,IRGC32399,IRGC418,IRGC55471,IRGC8191,IRGC38698,IRGC11010,IRGC17757,IRGC328,IRGC43325,IRGC43675,IRGC50448,IRGC66756,IRGC8244,IRGC26872,IRGC12793,IRGC38994,IRGC9060,IRGC9062,RA4952,IRGC31856' \
        --plan 'name=OsatJap_n24_OsatInd;t=OsatJap;qs=IRGC1107,IRGC2540,IRGC27630,IRGC32399,IRGC418,IRGC55471,IRGC8191,IRGC38698,IRGC11010,IRGC17757,IRGC328,IRGC43325,IRGC43675,IRGC50448,IRGC66756,IRGC8244,IRGC26872,IRGC12793,IRGC38994,IRGC9060,IRGC9062,RA4952,IRGC31856,OsatInd;o=OsatInd' \
        -y
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/rice82_test.yml

    sh 01_file.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_OsatJap_n8_OsatInd.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_OsatJap_n24_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_OsatJap_n24_OsatInd.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

# *Drosophila* Population Genomics Project (dpgp)

1. Create data.yml manually.

```bash
mkdir -p ~/data/alignment/dpgp82
cd ~/data/alignment/dpgp82

cat <<EOF > dpgp82_data.yml
---
data:
  - coverage: 37.8
    name: CK1
    original_id: 7227
  - coverage: 41.3
    name: CO15N
    original_id: 7227
  - coverage: 36.23
    name: ED10N
    original_id: 7227
  - coverage: 33.81
    name: EZ5N
    original_id: 7227
  - coverage: 38.44
    name: FR217
    original_id: 7227
  - coverage: 47.37
    name: GA185
    original_id: 7227
  - coverage: 37.08
    name: GU10
    original_id: 7227
  - coverage: 37.05
    name: KN6
    original_id: 7227
  - coverage: 33
    name: KR39
    original_id: 7227
  - coverage: 35.21
    name: KT1
    original_id: 7227
  - coverage: 36.97
    name: NG3N
    original_id: 7227
  - coverage: 28.41
    name: RC1
    original_id: 7227
  - coverage: 41.3
    name: RG15
    original_id: 7227
  - coverage: 40.78
    name: SP254
    original_id: 7227
  - coverage: 33.76
    name: TZ8
    original_id: 7227
  - coverage: 36.39
    name: UG7
    original_id: 7227
  - coverage: 42.83
    name: UM526
    original_id: 7227
  - coverage: 39.92
    name: ZI268
    original_id: 7227
  - coverage: 42.65
    name: ZL130
    original_id: 7227
  - coverage: 38.76
    name: ZO12
    original_id: 7227
  - coverage: 39.11
    name: ZS37
    original_id: 7227
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/dpgp82
    cd ~/data/alignment/dpgp82

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i dpgp82_data.yml \
        -o ~/Scripts/withncbi/pop/dpgp82_test.yml \
        -d ~/data/alignment/others/dpgp \
        -m name \
        -r '*.fa' \
        --opt group_name=dpgp82 \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/dpgp82' \
        --dd ~/data/alignment/Ensembl \
        --download 'name=Dmel;taxon=7227;sciname=Drosophila melanogaster' \
        --download 'name=Dsim;taxon=7240;sciname=Drosophila simulans' \
        --plan 'name=Dmel_n22_pop;t=Dmel;qs=CK1,CO15N,ED10N,EZ5N,FR217,GA185,GU10,KN6,KR39,KT1,NG3N,RC1,RG15,SP254,TZ8,UG7,UM526,ZI268,ZL130,ZO12,ZS37' \
        --plan 'name=Dmel_n22_Dsim;t=Dmel;qs=CK1,CO15N,ED10N,EZ5N,FR217,GA185,GU10,KN6,KR39,KT1,NG3N,RC1,RG15,SP254,TZ8,UG7,UM526,ZI268,ZL130,ZO12,ZS37,Dsim;o=Dsim' \
        -y
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/dpgp82_test.yml

    sh 01_file.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Dmel_n22_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Dmel_n22_Dsim.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

# Primates

1. Create an empty data.yml.

```bash
mkdir -p ~/data/alignment/primates82
cd ~/data/alignment/primates82

cat <<EOF > primates82_data.yml
---
data: []
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/primates82
    cd ~/data/alignment/primates82

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i primates82_data.yml \
        -o ~/Scripts/withncbi/pop/primates82_test.yml \
        --opt group_name=primates82 \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/primates82' \
        --opt phylo_tree='~/data/alignment/primates_13way.nwk' \
        --dd ~/data/alignment/Ensembl \
        --download 'name=Human;taxon=9606;sciname=Homo sapiens' \
        --download 'name=Chimp;taxon=9598;sciname=Pan troglodytes;coverage=6x sanger' \
        --download 'name=Gorilla;taxon=9595;sciname=Gorilla gorilla;coverage=2.1x sanger, 35x solexa' \
        --download 'name=Orangutan;taxon=9601;sciname=Pongo abelii;coverage=6x sanger' \
        --download 'name=Rhesus;taxon=9544;sciname=Macaca mulatta;coverage=6.1x sanger' \
        --plan 'name=HC_R;t=Human;qs=Chimp,Rhesus;o=Rhesus' \
        --plan 'name=HCGO_R;t=Human;qs=Chimp,Gorilla,Orangutan,Rhesus;o=Rhesus' \
        -y
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/primates82_test.yml

    sh 01_file.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_HC_R.sh
    sh 5_multi_cmd.sh

    sh plan_HCGO_R.sh
    sh 5_multi_cmd.sh
    ```

# Human individuals from Simons project

1. Create data.yml manually.

```bash
mkdir -p ~/data/alignment/human_simons
cd ~/data/alignment/human_simons

cat <<EOF > human_simons_data.yml
---
data:
  - name: HGDP00456
    population: Mbuti
    gender: M
    original_id: 9606
  - name: HGDP00521
    population: French
    gender: M
    original_id: 9606
  - name: HGDP00542
    population: Papuan
    gender: M
    original_id: 9606
  - name: HGDP00665
    population: Sardinian
    gender: M
    original_id: 9606
  - name: HGDP00778
    population: Han
    gender: M
    original_id: 9606
  - name: HGDP00927
    population: Yoruba
    gender: M
    original_id: 9606
  - name: HGDP00998
    population: Karitiana
    gender: M
    original_id: 9606
  - name: HGDP01029
    population: San
    gender: M
    original_id: 9606
  - name: HGDP01284
    population: Mandenka
    gender: M
    original_id: 9606
  - name: HGDP01307
    population: Dai
    gender: M
    original_id: 9606
  - name: MIXE
    population: Mixe
    gender: F
    original_id: 9606
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/human_simons
    cd ~/data/alignment/human_simons

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i human_simons_data.yml \
        -o ~/Scripts/withncbi/pop/human_simons_test.yml \
        -d ~/data/alignment/others/simons-phase0 \
        -m name \
        -r '*.fa.gz' \
        --opt group_name=human_simons \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/human_simons' \
        --dd ~/data/alignment/Ensembl \
        --download 'name=Human;taxon=9606;sciname=Homo sapiens' \
        --download 'name=Chimp;taxon=9598;sciname=Pan troglodytes;coverage=6x sanger' \
        --plan 'name=Human_n12_pop;t=Human;qs=HGDP00456,HGDP00521,HGDP00542,HGDP00665,HGDP00778,HGDP00927,HGDP00998,HGDP01029,HGDP01284,HGDP01307,MIXE' \
        --plan 'name=Human_n12_Chimp;t=Human;qs=HGDP00456,HGDP00521,HGDP00542,HGDP00665,HGDP00778,HGDP00927,HGDP00998,HGDP01029,HGDP01284,HGDP01307,MIXE,Chimp;o=Chimp' \
        -y
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    cd ~/data/alignment/human_simons
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 23 -i ~/Scripts/withncbi/pop/human_simons_test.yml

    sh 01_file.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Human_n12_pop.sh
    sh 5_multi_cmd.sh

    sh plan_Human_n12_Chimp.sh
    sh 5_multi_cmd.sh
    ```

# *Caenorhabditis elegans* million mutation project (cele_mmp)

1. Create data.yml manually.

    Information came from
    [here](http://genome.cshlp.org/content/suppl/2013/08/20/gr.157651.113.DC2/Supplemental_Table_3.xls).

```bash
mkdir -p ~/data/alignment/cele82
cd ~/data/alignment/cele82

cat <<EOF > cele82_data.yml
---
data:
  - coverage: 29
    name: AB1
    original_id: 6239
  - coverage: 27.6
    name: AB3
    original_id: 6239
  - coverage: 28.8
    name: CB4853
    original_id: 6239
  - coverage: 29.9
    name: CB4854
    original_id: 6239
  - coverage: 30.8
    name: CB4856
    original_id: 6239
  - coverage: 30.1
    name: ED3017
    original_id: 6239
  - coverage: 29.2
    name: ED3021
    original_id: 6239
  - coverage: 28.8
    name: ED3040
    original_id: 6239
  - coverage: 29.2
    name: ED3042
    original_id: 6239
  - coverage: 28.6
    name: ED3049
    original_id: 6239
  - coverage: 28.6
    name: ED3052
    original_id: 6239
  - coverage: 30.3
    name: ED3057
    original_id: 6239
  - coverage: 28.9
    name: ED3072
    original_id: 6239
  - coverage: 29.4
    name: GXW1
    original_id: 6239
  - coverage: 30.4
    name: JU1088
    original_id: 6239
  - coverage: 28.1
    name: JU1171
    original_id: 6239
  - coverage: 29.5
    name: JU1400
    original_id: 6239
  - coverage: 28.9
    name: JU1401
    original_id: 6239
  - coverage: 27.9
    name: JU1652
    original_id: 6239
  - coverage: 28.5
    name: JU258
    original_id: 6239
  - coverage: 28.7
    name: JU263
    original_id: 6239
  - coverage: 30
    name: JU300
    original_id: 6239
  - coverage: 30.6
    name: JU312
    original_id: 6239
  - coverage: 28.2
    name: JU322
    original_id: 6239
  - coverage: 30.3
    name: JU345
    original_id: 6239
  - coverage: 29.2
    name: JU360
    original_id: 6239
  - coverage: 29.3
    name: JU361
    original_id: 6239
  - coverage: 30.9
    name: JU394
    original_id: 6239
  - coverage: 27.8
    name: JU397
    original_id: 6239
  - coverage: 28.6
    name: JU533
    original_id: 6239
  - coverage: 28.5
    name: JU642
    original_id: 6239
  - coverage: 30.2
    name: JU775
    original_id: 6239
  - coverage: 30.1
    name: KR314
    original_id: 6239
  - coverage: 29.5
    name: LKC34
    original_id: 6239
  - coverage: 28.7
    name: MY1
    original_id: 6239
  - coverage: 25.9
    name: MY14
    original_id: 6239
  - coverage: 24
    name: MY16
    original_id: 6239
  - coverage: 29
    name: MY2
    original_id: 6239
  - coverage: 24.1
    name: MY6
    original_id: 6239
  - coverage: 28.8
    name: PX174
    original_id: 6239
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/cele82
    cd ~/data/alignment/cele82

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i cele82_data.yml \
        -o ~/Scripts/withncbi/pop/cele82_test.yml \
        -d ~/data/alignment/others/cele \
        -m name \
        -r '*.vcf.fasta' \
        --opt group_name=cele82 \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/cele82' \
        --opt rm_species='Caenorhabditis elegans' \
        --dd ~/data/alignment/Ensembl \
        --download 'name=Cele;taxon=6239;sciname=Caenorhabditis elegans' \
        --plan 'name=Cele_n41_pop;t=Cele;qs=AB1,AB3,CB4853,CB4854,CB4856,ED3017,ED3021,ED3040,ED3042,ED3049,ED3052,ED3057,ED3072,GXW1,JU1088,JU1171,JU1400,JU1401,JU1652,JU258,JU263,JU300,JU312,JU322,JU345,JU360,JU361,JU394,JU397,JU533,JU642,JU775,KR314,LKC34,MY14,MY16,MY1,MY2,MY6,PX174' \
        -y
    ```

2. Rest routing things.

    ```bash
    cd ~/data/alignment/cele82

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/cele82_test.yml

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

# *Dictyostelium* WGS

1. RM species

    It't OK to not specifies RM species. Protists have very few repeats records.

    ```bash
    /usr/local/Cellar/repeatmasker/4.0.5/libexec/util/queryTaxonomyDatabase.pl -taxDBFile /usr/local/Cellar/repeatmasker/4.0.5/libexec/Libraries/taxonomy.dat -species Dictyosteliida
    /usr/local/Cellar/repeatmasker/4.0.5/libexec/util/queryRepeatDatabase.pl -species Dictyosteliida -stat
    ```

2. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/Protists/dictyostelium
    cd ~/data/alignment/Protists/dictyostelium

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Protists/GENOMES/dictyostelium/WGS/dictyostelium.data.yml \
        -o ~/Scripts/withncbi/pop/dictyostelium_test.yml \
        -d ~/data/alignment/Protists/GENOMES/dictyostelium/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=dictyostelium \
        --opt base_dir='~/data/alignment/Protists' \
        --opt data_dir='~/data/alignment/Protists/dictyostelium' \
        --dd ~/data/alignment/Protists/GENOMES/dictyostelium/DOWNLOAD \
        --download 'name=Ddis_AX4;taxon=352472;sciname=Dictyostelium discoideum AX4' \
        --download 'name=Dictyostelium_purpureum;taxon=5786;sciname=Dictyostelium purpureum' \
        --download 'name=Dictyostelium_firmibasis;taxon=79012;sciname=Dictyostelium firmibasis' \
        --download 'name=Dictyostelium_fasciculatum;taxon=261658;sciname=Dictyostelium fasciculatum' \
        --download 'name=Dictyostelium_citrinum;taxon=361072;sciname=Dictyostelium citrinum' \
        --download 'name=Dictyostelium_intermedium;taxon=361076;sciname=Dictyostelium intermedium' \
        -y
    ```

3. Rest routing things.

    ```bash
    cd ~/data/alignment/Protists/dictyostelium

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i ~/Scripts/withncbi/pop/dictyostelium_test.yml

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

4. Pick outgroup.

    Dictyostelium_citrinum or Dictyostelium_firmibasis assemblies.

# *Dictyostelium discoideum*

1. Create data.yml manually.

    Coverages are not real.

```bash
mkdir -p ~/data/alignment/Protists/Ddis
cd ~/data/alignment/Protists/Ddis

cat <<EOF > ddis_data.yml
---
data:
  - coverage: 10
    name: 68
    original_id: 44689
  - coverage: 10
    name: 70
    original_id: 44689
  - coverage: 10
    name: QS11
    original_id: 44689
  - coverage: 10
    name: QS17
    original_id: 44689
  - coverage: 10
    name: QS18
    original_id: 44689
  - coverage: 10
    name: QS23
    original_id: 44689
  - coverage: 10
    name: QS36
    original_id: 44689
  - coverage: 10
    name: QS37
    original_id: 44689
  - coverage: 10
    name: QS4
    original_id: 44689
  - coverage: 10
    name: QS69
    original_id: 44689
  - coverage: 10
    name: QS73
    original_id: 44689
  - coverage: 10
    name: QS74
    original_id: 44689
  - coverage: 10
    name: QS80
    original_id: 44689
  - coverage: 10
    name: QS9
    original_id: 44689
  - coverage: 10
    name: S224
    original_id: 44689
  - coverage: 10
    name: WS14
    original_id: 44689
  - coverage: 10
    name: WS15
    original_id: 44689
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/Protists/Ddis
    cd ~/data/alignment/Protists/Ddis

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ddis_data.yml \
        -o ~/Scripts/withncbi/pop/ddis_test.yml \
        -d ~/data/alignment/others/dicty \
        -m name \
        -r '*.fa' \
        --opt group_name=Ddis \
        --opt base_dir='~/data/alignment/Protists' \
        --opt data_dir='~/data/alignment/Protists/Ddis' \
        --dd ~/data/alignment/Protists/GENOMES/Ddis/DOWNLOAD \
        --download 'name=AX4;taxon=352472;sciname=Dictyostelium discoideum AX4' \
        --download 'name=Dfir;taxon=79012;sciname=Dictyostelium firmibasis' \
        --download 'name=Dcit;taxon=361072;sciname=Dictyostelium citrinum' \
        --plan 'name=Ddis_n18_pop;t=AX4;qs=68,70,QS11,QS17,QS18,QS23,QS36,QS37,QS4,QS69,QS73,QS74,QS80,QS9,S224,WS14,WS15' \
        --plan 'name=Ddis_n11_pop;t=AX4;qs=68,70,QS23,QS36,QS37,QS69,QS73,QS74,QS80,S224' \
        --plan 'name=Ddis_n11_Dfir;t=AX4;qs=68,70,QS23,QS36,QS37,QS69,QS73,QS74,QS80,S224,Dfir;o=Dfir' \
        --plan 'name=Ddis_n11_Dcit;t=AX4;qs=68,70,QS23,QS36,QS37,QS69,QS73,QS74,QS80,S224,Dcit;o=Dcit' \
        -y
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/ddis_test.yml

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
    sh plan_Ddis_n11_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Ddis_n18_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Ddis_n11_Dfir.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Ddis_n11_Dcit.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```


# Mouse

1. Create data.yml manually.

```bash
mkdir -p ~/data/alignment/mouse82
cd ~/data/alignment/mouse82

cat <<EOF > mouse82_data.yml
---
data:
  - name: 129S1_SvImJ
    original_id: 10090
  - name: A_J
    original_id: 10090
  - name: AKR_J
    original_id: 10090
  - name: BALB_cJ
    original_id: 10090
  - name: C3H_HeJ
    original_id: 10090
  - name: C57BL_6NJ
    original_id: 10090
  - name: CAROLI_EiJ
    original_id: 10090
  - name: CAST_EiJ
    taxon: 10091
    sciname: Mus musculus castaneus
  - name: CBA_J
    original_id: 10090
  - name: DBA_2J
    original_id: 10090
  - name: FVB_NJ
    original_id: 10090
  - name: LP_J
    original_id: 10090
  - name: NOD_ShiLtJ
    original_id: 10090
  - name: NZO_HlLtJ
    original_id: 10090
  - name: Pahari_EiJ
    original_id: 10090
  - name: PWK_PhJ
    taxon: 39442
    sciname: Mus musculus musculus
  - name: SPRET_EiJ
    taxon: 10096
    sciname: Mus spretus
  - name: WSB_EiJ
    taxon: 10092
    sciname: Mus musculus domesticus
EOF
```

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/mouse82
    cd ~/data/alignment/mouse82

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i mouse82_data.yml \
        -o ~/Scripts/withncbi/pop/mouse82_test.yml \
        -d ~/data/alignment/others/sanger-mouse \
        -m name \
        -r '*.fa.masked.gz' \
        --opt group_name=mouse82 \
        --opt base_dir='~/data/alignment' \
        --opt data_dir='~/data/alignment/mouse82' \
        --dd ~/data/alignment/Ensembl \
        --download 'name=Mouse;taxon=10090;sciname=Mus musculus' \
        --download 'name=Rat;taxon=10116;sciname=Rattus norvegicus' \
        --plan 'name=Mouse_n11_pop;t=Mouse;qs=129S1_SvImJ,A_J,AKR_J,BALB_cJ,C3H_HeJ,CBA_J,DBA_2J,LP_J,NOD_ShiLtJ,NZO_HlLtJ' \
        --plan 'name=Mouse_n11_SPRET_EiJ;t=Mouse;qs=129S1_SvImJ,A_J,AKR_J,BALB_cJ,C3H_HeJ,CBA_J,DBA_2J,LP_J,NOD_ShiLtJ,NZO_HlLtJ,SPRET_EiJ;o=SPRET_EiJ' \
        --plan 'name=Mouse_n14_pop;t=Mouse;qs=129S1_SvImJ,A_J,AKR_J,BALB_cJ,C3H_HeJ,CAROLI_EiJ,CBA_J,DBA_2J,FVB_NJ,LP_J,NOD_ShiLtJ,NZO_HlLtJ,Pahari_EiJ' \
        --plan 'name=Mouse_n14_SPRET_EiJ;t=Mouse;qs=129S1_SvImJ,A_J,AKR_J,BALB_cJ,C3H_HeJ,CAROLI_EiJ,CBA_J,DBA_2J,FVB_NJ,LP_J,NOD_ShiLtJ,NZO_HlLtJ,Pahari_EiJ,SPRET_EiJ;o=SPRET_EiJ' \
        --plan 'name=Mouse_n5_Rat;t=Mouse;qs=CAST_EiJ,PWK_PhJ,SPRET_EiJ,WSB_EiJ,Rat;o=Rat' \
        -y
    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i ~/Scripts/withncbi/pop/mouse82_test.yml

    sh 01_file.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Mouse_n11_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Mouse_n11_SPRET_EiJ.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Mouse_n5_Rat.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh    
    ```

