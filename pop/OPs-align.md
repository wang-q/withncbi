# Aligning steps for each groups

Less detailed than Trichoderma in [README.md](README.md), but include examples for genomes out of
WGS, which usually in better assembling levels.

[TOC levels=1-3]: # ""

- [Aligning steps for each groups](#aligning-steps-for-each-groups)
- [*Arabidopsis* 19 genomes](#arabidopsis-19-genomes)
- [*Orazy sativa* Japonica 24 genomes](#orazy-sativa-japonica-24-genomes)
- [*Drosophila* Population Genomics Project (dpgp)](#drosophila-population-genomics-project-dpgp)
- [Primates](#primates)
- [Human individuals from Simons project](#human-individuals-from-simons-project)
- [*Caenorhabditis elegans* million mutation project (cele_mmp)](#caenorhabditis-elegans-million-mutation-project-cele_mmp)
- [*Dictyostelium* WGS](#dictyostelium-wgs)
- [*Dictyostelium discoideum*](#dictyostelium-discoideum)
- [Mouse](#mouse)


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

