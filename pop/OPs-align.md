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

### *Scer_new* WGS and ASSEMBLY

1. `gen_pop_conf.pl`

    ```bash
	# create downloaded list
    cat ~/data/alignment/Fungi/GENOMES/scer_new/DOWNLOAD/scer_new.seq.csv \
        | grep -v "^#" \
        | cut -d',' -f1,3 \
        | uniq \
        | perl -nl -a -F"," -e 'printf qq{    --download "name=%s;taxon=%s" \\\n}, $F[0], $F[1];'

    mkdir -p ~/data/alignment/Fungi/scer_new
    cd ~/data/alignment/Fungi/scer_new

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/Fungi/GENOMES/scer_new/WGS/scer_new.data.yml \
        -o ~/Scripts/withncbi/pop/scer_new_test.yml \
        -d ~/data/alignment/Fungi/GENOMES/scer_new/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=scer_new \
        --opt base_dir='~/data/alignment/Fungi' \
        --opt data_dir='~/data/alignment/Fungi/scer_new' \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/Fungi/GENOMES/scer_new/DOWNLOAD \
	    --download "name=S288c;taxon=559292" \
	    --download "name=EC1118;taxon=643680" \
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
        --plan 'name=three_way;t=S288c;qs=RM11_1a,YJM789'

    ```

2. Rest routing things.

    ```bash
    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 8 -i ~/Scripts/withncbi/pop/scer_new_test.yml

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
    sh plan_three_way.sh

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
