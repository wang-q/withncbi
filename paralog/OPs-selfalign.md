# Self-aligning steps for each groups

## Taxonomy for Ensembl species

```bash
mkdir -p ~/data/alignment/self
cd ~/data/alignment/self

perl        ~/Scripts/withncbi/taxon/strain_info.pl \
    --id    559292                                  --name 559292=S288c  \
    --id    7227                                    --name 7227=Dmel     \
    --id    9606                                    --name 9606=Human    \
    --id    6239                                    --name 6239=Cele     \
    --id    352472                                  --name 352472=Ddis   \
    --id    10090                                   --name 10090=Mouse   \
    --id    3880                                    --name 3880=Mtru     \
    --id    3847                                    --name 3847=Gmax     \
    --id    3712                                    --name 3712=Bole     \
    --id    3711                                    --name 3711=Brap     \
    --id    3702                                    --name 3702=Atha     \
    --id    59689                                   --name 59689=Alyr    \
    --id    29760                                   --name 29760=Vvin    \
    --id    4081                                    --name 4081=Slyc     \
    --id    4113                                    --name 4113=Stub     \
    --id    4641                                    --name 4641=Macu     \
    --id    4555                                    --name 4555=Sita     \
    --id    4558                                    --name 4558=Sbic     \
    --id    39947                                   --name 39947=OsatJap \
    --id    15368                                   --name 15368=Bdis    \
    --file ensembl_taxon.csv                        \
    --entrez

```

### Plants

| Name                       | Classification | Taxon ID | Used |
| :---                       | :---           | :---     | :--- |
| Aegilops tauschii          | Liliopsida     | 37682    |      |
| Amborella trichopoda       | Amborellales   | 13333    | -    |
| Arabidopsis lyrata         | eudicotyledons | 81972    | -    |
| Arabidopsis thaliana       | eudicotyledons | 3702     | -    |
| Brachypodium distachyon    | Liliopsida     | 15368    | -    |
| Brassica oleracea          | eudicotyledons | 109376   |      |
| Brassica rapa              | eudicotyledons | 51351    | -    |
| Chlamydomonas reinhardtii  | Chlorophyta    | 3055     |      |
| Cyanidioschyzon merolae    | Rhodophyta     | 280699   |      |
| Glycine max                | eudicotyledons | 3847     | -    |
| Hordeum vulgare            | Liliopsida     | 112509   |      |
| Leersia perrieri           | Liliopsida     | 77586    |      |
| Medicago truncatula        | eudicotyledons | 3880     | -    |
| Musa acuminata             | Liliopsida     | 214687   | -    |
| Oryza barthii              | Liliopsida     | 65489    |      |
| Oryza brachyantha          | Liliopsida     | 4533     |      |
| Oryza glaberrima           | Liliopsida     | 4538     |      |
| Oryza glumaepatula         | Liliopsida     | 40148    |      |
| Oryza longistaminata       | Liliopsida     | 4528     |      |
| Oryza meridionalis         | Liliopsida     | 40149    |      |
| Oryza nivara               | Liliopsida     | 4536     |      |
| Oryza punctata             | Liliopsida     | 4537     |      |
| Oryza rufipogon            | Liliopsida     | 4529     |      |
| Oryza sativa Indica        | Liliopsida     | 39946    |      |
| Oryza sativa Japonica      | Liliopsida     | 39947    | -    |
| Ostreococcus lucimarinus   | Chlorophyta    | 436017   |      |
| Physcomitrella patens      | Bryophyta      | 3218     |      |
| Populus trichocarpa        | eudicotyledons | 3694     |      |
| Prunus persica             | eudicotyledons | 3760     |      |
| Selaginella moellendorffii | Lycopodiophyta | 88036    |      |
| Setaria italica            | Liliopsida     | 4555     | -    |
| Solanum lycopersicum       | eudicotyledons | 4081     | -    |
| Solanum tuberosum          | eudicotyledons | 4113     | -    |
| Sorghum bicolor            | Liliopsida     | 4558     | -    |
| Theobroma cacao            | eudicotyledons | 3641     |      |
| Triticum aestivum          | Liliopsida     | 4565     |      |
| Triticum urartu            | Liliopsida     | 4572     |      |
| Vitis vinifera             | eudicotyledons | 29760    | -    |
| Zea mays                   | Liliopsida     | 4577     |      |

## Yeast S288c

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name yeast \
    --parallel 8 \
    -t S288c

cd ~/data/alignment/self/yeast

bash 1_real_chr.sh
bash 3_self_cmd.sh
time bash 4_proc_cmd.sh
# real    0m46.075s
# user    1m19.171s
# sys     0m42.335s
bash 5_circos_cmd.sh
```

## Arabidopsis

### 1. full chromosomes

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name arabidopsis \
    --parallel 8 \
    -t Atha

cd ~/data/alignment/self/arabidopsis

# ath centromere position from the follow file:
# ftp://ftp.arabidopsis.org/home/tair/Sequences/whole_chromosomes/tair9_Assembly_gaps.gff
# chr1 14511722	14803970
# chr2 3611839	3633423
# chr3 13589757	13867121
# chr4 3133664	3133674
# chr5 11194538	11723210

# Hosouchi, T., Kumekawa, N., Tsuruoka, H. & Kotani, H. Physical Map-Based Sizes of the Centromeric Regions of Arabidopsis thaliana Chromosomes 1, 2, and 3. DNA Res 9, 117®C121 (2002).
TAB=$'\t'
cat <<EOF > Atha.kary.tsv
#chrom${TAB}chromStart${TAB}chromEnd${TAB}name${TAB}gieStain
1${TAB}1${TAB}14200000${TAB}p1${TAB}gpos50${TAB}#14.2M
1${TAB}14200000${TAB}15627671${TAB}p1${TAB}acen
1${TAB}15627671${TAB}30427671${TAB}q1${TAB}gpos50${TAB}#14.8M
2${TAB}1${TAB}3000000${TAB}p2${TAB}gpos50${TAB}#3.0M
2${TAB}3000000${TAB}3898289${TAB}p2${TAB}acen
2${TAB}3898289${TAB}19698289${TAB}q2${TAB}gpos50${TAB}#15.8M
3${TAB}1${TAB}13200000${TAB}p3${TAB}gpos50${TAB}#13.2M
3${TAB}13200000${TAB}14459830${TAB}p3${TAB}acen
3${TAB}14459830${TAB}23459830${TAB}q3${TAB}gpos50${TAB}#9M
4${TAB}1${TAB}3000000${TAB}p4${TAB}gpos50${TAB}#3M
4${TAB}3000000${TAB}5085056${TAB}p4${TAB}acen
4${TAB}5085056${TAB}18585056${TAB}q4${TAB}gpos50${TAB}#13.5M
5${TAB}1${TAB}11100000${TAB}p5${TAB}gpos50${TAB}#${TAB}1.1M
5${TAB}11100000${TAB}12575502${TAB}p5${TAB}acen
5${TAB}12575502${TAB}26975502${TAB}q5${TAB}gpos50${TAB}#14.4M
EOF
bash ~/share/circos/data/karyotype/parse.karyotype Atha.kary.tsv > Processing/Atha/karyotype.Atha.txt

bash 1_real_chr.sh
time bash 3_self_cmd.sh
# real    25m15.804s
# user    161m4.543s
# sys     1m10.534s
time bash 4_proc_cmd.sh
# real    10m21.329s
# user    27m48.637s
# sys     12m24.249s
bash 5_circos_cmd.sh
```

### 2. partition sequences

```bash
cd ~/data/alignment/self

mkdir -p ~/data/alignment/self/arabidopsis_parted/Genomes
perl ~/Scripts/egaz/part_seq.pl \
    -i ~/data/alignment/Ensembl/Atha \
    -o ~/data/alignment/self/arabidopsis_parted/Genomes/Atha \
    --chunk 10010000 --overlap 10000

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --parted \
    --length 1000  \
    --use_name \
    --norm \
    --name arabidopsis_parted \
    --parallel 8 \
    -t Atha

cd ~/data/alignment/self/arabidopsis_parted

bash 1_real_chr.sh
time bash 3_self_cmd.sh
# real    21m10.875s
# user    156m48.601s
# sys     1m54.846s
time bash 4_proc_cmd.sh
# real    9m33.086s
# user    24m56.361s
# sys     11m21.288s
bash 5_circos_cmd.sh
```

### 3. Comparison

**11.50% vs 9.96%. Use full chromosomes if the running time is acceptable.**

## Cele

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name worm \
    --parallel 8 \
    -t Cele

cd ~/data/alignment/self/worm

bash 1_real_chr.sh
time bash 3_self_cmd.sh
# real    26m9.090s
# user    150m18.730s
# sys     0m49.323s
time bash 4_proc_cmd.sh
# real    3m24.714s
# user    6m42.150s
# sys     4m1.954s
bash 5_circos_cmd.sh
```

## Ddis

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name dicty \
    --parallel 8 \
    -t Ddis

cd ~/data/alignment/self/dicty

bash 1_real_chr.sh
time bash 3_self_cmd.sh
# real    1m53.391s
# user    13m25.636s
# sys     0m7.271s
time bash 4_proc_cmd.sh
# real    353m10.864s
# user    364m44.545s
# sys     2m29.146s
bash 5_circos_cmd.sh
```

## Human

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name human \
    --parallel 8 \
    -t Human

cd ~/data/alignment/self/human

bash 1_real_chr.sh
time bash 3_self_cmd.sh
# real    4156m4.259s
# user    32775m39.086s
# sys     227m18.093s
time bash 4_proc_cmd.sh
# real    940m13.950s
# user    3947m24.734s
# sys     248m12.852s
bash 5_circos_cmd.sh
```

## Mouse

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name mouse \
    --parallel 8 \
    -t Mouse

cd ~/data/alignment/self/mouse

bash 1_real_chr.sh
time bash 3_self_cmd.sh
# real    3801m39.031s
# user    29792m18.127s
# sys     133m41.809s
time bash 4_proc_cmd.sh
# real    1750m33.958s
# user    5020m19.404s
# sys     365m30.568s
bash 5_circos_cmd.sh
```

## All plants

### Full chromosomes

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name plants \
    --parallel 12 \
    -q Alyr \
    -q OsatJap \
    -q Sbic \
    -t Atha

cd ~/data/alignment/self/plants

bash 1_real_chr.sh
time bash 3_self_cmd.sh
# real    658m10.559s
# user    5623m25.407s
# sys     74m46.390s
time bash 4_proc_cmd.sh
# real    971m49.779s
# user    2130m57.202s
# sys     435m29.950s
bash 5_circos_cmd.sh
```

### Partitioned chromosomes

`Bole` contains exact matched pieces with copy number large than one thousand.

It turns out that poor assemblies tend to have this phenomenon.

To speed up processing, use partitioned sequences in `3_self_cmd.sh` and `--discard 50` in `4_proc_cmd.sh`.

```bash
cd ~/data/alignment/self

mkdir -p ~/data/alignment/self/plants_parted/Genomes

for name in Mtru Gmax Bole Brap Alyr Vvin Slyc Stub Macu Sita OsatJap Bdis Atha
do
    echo "==> ${name}"
    perl ~/Scripts/egaz/part_seq.pl \
        -i ~/data/alignment/Ensembl/${name} \
        -o ~/data/alignment/self/plants_parted/Genomes/${name} \
        --chunk 10010000 --overlap 10000
done

perl ~/Scripts/egaz/self_batch.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --length 1000  \
    --use_name \
    --norm \
    --name plants_parted \
    --parallel 12 \
    -q Mtru \
    -q Gmax \
    -q Bole \
    -q Brap \
    -q Alyr \
    -q Vvin \
    -q Slyc \
    -q Stub \
    -q OsatJap \
    -q Bdis \
    -t Atha \
    --parted

    # cost days
    # -q Macu \
    # -q Sita \

cd ~/data/alignment/self/plants_parted

# perl ~/Scripts/withncbi/ensembl/chr_kary.pl -e oryza_sativa_core_29_82_7
# bash ~/share/circos/data/karyotype/parse.karyotype oryza_sativa_core_29_82_7.kary.tsv > Processing/OsatJap/karyotype.OsatJap.txt

# ensembldb.ensembl.org         5306
# mysql-eg-publicsql.ebi.ac.uk  4157
# mysql -hmysql-eg-publicsql.ebi.ac.uk -P4157 -uanonymous
# perl ~/Scripts/withncbi/ensembl/chr_kary.pl -s mysql-eg-publicsql.ebi.ac.uk --port 4157 -u anonymous -p '' -e oryza_sativa_core_29_82_7

bash 1_real_chr.sh
bash 3_self_cmd.sh
bash 4_proc_cmd.sh
bash 5_circos_cmd.sh
```
