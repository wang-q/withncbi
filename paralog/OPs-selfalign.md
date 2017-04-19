# Self-aligning steps for each groups


[TOC levels=1-3]: # " "
- [Self-aligning steps for each groups](#self-aligning-steps-for-each-groups)
- [Taxonomy for Ensembl species](#taxonomy-for-ensembl-species)
    - [Plants](#plants)
- [Model organisms](#model-organisms)
    - [Ecoli K-12 MG1655](#ecoli-k-12-mg1655)
    - [Yeast S288c](#yeast-s288c)
    - [Dmel](#dmel)
    - [Cele](#cele)
    - [Ddis](#ddis)
    - [Human](#human)
    - [Mouse](#mouse)
    - [Arabidopsis](#arabidopsis)
        - [Atha: full chromosomes](#atha-full-chromosomes)
        - [Atha: partition sequences](#atha-partition-sequences)
- [All plants](#all-plants)
    - [Plants: full chromosomes](#plants-full-chromosomes)
    - [Plants: partitioned chromosomes](#plants-partitioned-chromosomes)


# Taxonomy for Ensembl species

```bash
mkdir -p ~/data/alignment/self
cd ~/data/alignment/self

perl        ~/Scripts/withncbi/taxon/strain_info.pl \
    --id    511145 --name 511145=MG1655 \
    --id    559292 --name 559292=S288c  \
    --id    7227   --name 7227=Dmel     \
    --id    9606   --name 9606=Human    \
    --id    6239   --name 6239=Cele     \
    --id    352472 --name 352472=Ddis   \
    --id    10090  --name 10090=Mouse   \
    --id    3880   --name 3880=Mtru     \
    --id    3847   --name 3847=Gmax     \
    --id    3712   --name 3712=Bole     \
    --id    3711   --name 3711=Brap     \
    --id    3702   --name 3702=Atha     \
    --id    59689  --name 59689=Alyr    \
    --id    29760  --name 29760=Vvin    \
    --id    4081   --name 4081=Slyc     \
    --id    4113   --name 4113=Stub     \
    --id    4641   --name 4641=Macu     \
    --id    4555   --name 4555=Sita     \
    --id    4558   --name 4558=Sbic     \
    --id    39947  --name 39947=OsatJap \
    --id    15368  --name 15368=Bdis    \
    --file ensembl_taxon.csv            \
    --entrez
```

## Plants

| Name                       | Classification | Taxon ID | Used |
|:---------------------------|:---------------|:---------|:----:|
| Aegilops tauschii          | Liliopsida     | 37682    |      |
| Amborella trichopoda       | Amborellales   | 13333    |  -   |
| Arabidopsis lyrata         | eudicotyledons | 81972    |  -   |
| Arabidopsis thaliana       | eudicotyledons | 3702     |  -   |
| Brachypodium distachyon    | Liliopsida     | 15368    |  -   |
| Brassica oleracea          | eudicotyledons | 109376   |      |
| Brassica rapa              | eudicotyledons | 51351    |  -   |
| Chlamydomonas reinhardtii  | Chlorophyta    | 3055     |      |
| Cyanidioschyzon merolae    | Rhodophyta     | 280699   |      |
| Glycine max                | eudicotyledons | 3847     |  -   |
| Hordeum vulgare            | Liliopsida     | 112509   |      |
| Leersia perrieri           | Liliopsida     | 77586    |      |
| Medicago truncatula        | eudicotyledons | 3880     |  -   |
| Musa acuminata             | Liliopsida     | 214687   |  -   |
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
| Oryza sativa Japonica      | Liliopsida     | 39947    |  -   |
| Ostreococcus lucimarinus   | Chlorophyta    | 436017   |      |
| Physcomitrella patens      | Bryophyta      | 3218     |      |
| Populus trichocarpa        | eudicotyledons | 3694     |      |
| Prunus persica             | eudicotyledons | 3760     |      |
| Selaginella moellendorffii | Lycopodiophyta | 88036    |      |
| Setaria italica            | Liliopsida     | 4555     |  -   |
| Solanum lycopersicum       | eudicotyledons | 4081     |  -   |
| Solanum tuberosum          | eudicotyledons | 4113     |  -   |
| Sorghum bicolor            | Liliopsida     | 4558     |  -   |
| Theobroma cacao            | eudicotyledons | 3641     |      |
| Triticum aestivum          | Liliopsida     | 4565     |      |
| Triticum urartu            | Liliopsida     | 4572     |      |
| Vitis vinifera             | eudicotyledons | 29760    |  -   |
| Zea mays                   | Liliopsida     | 4577     |      |

TODO:

* Sequencing strategy and coverage: BAC, WGS, Illumina
    * N50 of contigs and scaffolds

# Model organisms

## Ecoli K-12 MG1655

Pretend to be an Ensembl project.

Use GenBank accessions instead of RefSeq accessions.

```bash
mkdir -p ~/data/alignment/Ensembl/MG1655
cd ~/data/alignment/Ensembl/MG1655

#curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_000913.3&rettype=fasta&retmode=txt" \
#    > NC_000913.fa

cat <<EOF > name_seq.csv
strain_name,accession,strain_taxon_id,seq_name # This line is needed
MG1655,U00096,
EOF

perl ~/Scripts/withncbi/taxon/batch_get_seq.pl -p -f name_seq.csv

mv MG1655/* .
rm -fr MG1655

# simplify header
perl -pi -e 's/^>(\w+).+/>$1/' U00096.fasta

RepeatMasker U00096.fasta -xsmall -gff -parallel 8

mv U00096.fasta.masked U00096.fa
mv U00096.fasta.out.gff U00096.rm.gff
rm U00096.fasta*

# create anno.yml
perl ~/Scripts/alignDB/util/gff2anno.pl \
    --type CDS --remove \
    U00096.gff \
    > cds.yml

perl ~/Scripts/alignDB/util/gff2anno.pl \
    --remove \
    U00096.rm.gff \
    > repeat.yml

runlist merge \
    repeat.yml \
    cds.yml \
    -o anno.yml

rm repeat.yml cds.yml

```

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name ecoli \
    --parallel 8 \
    -t MG1655

bash ecoli/1_real_chr.sh
bash ecoli/3_self_cmd.sh
bash ecoli/4_proc_cmd.sh
bash ecoli/5_circos_cmd.sh
```

## Yeast S288c

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name yeast \
    --parallel 8 \
    -t S288c

bash yeast/1_real_chr.sh
bash yeast/3_self_cmd.sh
# real    0m46.075s
time bash yeast/4_proc_cmd.sh
bash yeast/5_circos_cmd.sh
```

## Dmel

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name fly \
    --parallel 8 \
    -t Dmel

bash fly/1_real_chr.sh
#real    18m1.867s
#user    126m31.346s
#sys     0m21.494s
time bash fly/3_self_cmd.sh
#real    8m45.703s
#user    18m48.943s
#sys     7m3.803s
time bash fly/4_proc_cmd.sh
bash fly/5_circos_cmd.sh
```

## Cele

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name worm \
    --parallel 8 \
    -t Cele

bash worm/1_real_chr.sh
#real    16m41.346s
#user    116m6.311s
#sys     0m16.111s
time bash worm/3_self_cmd.sh
#real    4m56.288s
#user    9m59.978s
#sys     5m41.933s
time bash worm/4_proc_cmd.sh
bash worm/5_circos_cmd.sh
```

## Ddis

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name dicty \
    --parallel 8 \
    -t Ddis

bash dicty/1_real_chr.sh
# real    1m53.391s
time bash dicty/3_self_cmd.sh
# real    353m10.864s
time bash dicty/4_proc_cmd.sh
bash dicty/5_circos_cmd.sh
```

## Human

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name human \
    --parallel 8 \
    -t Human

cd ~/data/alignment/self/human

bash 1_real_chr.sh
# real    4156m4.259s
time bash 3_self_cmd.sh
# real    940m13.950s
time bash 4_proc_cmd.sh
bash 5_circos_cmd.sh
```

## Mouse

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name mouse \
    --parallel 8 \
    -t Mouse

cd ~/data/alignment/self/mouse

bash 1_real_chr.sh
# real    3801m39.031s
time bash 3_self_cmd.sh
# real    1750m33.958s
time bash 4_proc_cmd.sh
bash 5_circos_cmd.sh
```

## Arabidopsis

Comparison: **11.50% vs 9.96%**. Use full chromosomes if the running time is acceptable.

### Atha: full chromosomes

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
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
5${TAB}1${TAB}11100000${TAB}p5${TAB}gpos50${TAB}#1.1M
5${TAB}11100000${TAB}12575502${TAB}p5${TAB}acen
5${TAB}12575502${TAB}26975502${TAB}q5${TAB}gpos50${TAB}#14.4M
EOF

mkdir -p Processing/Atha
bash ~/share/circos/data/karyotype/parse.karyotype Atha.kary.tsv > Processing/Atha/karyotype.Atha.txt

bash 1_real_chr.sh
# real    25m15.804s
time bash 3_self_cmd.sh
# real    10m21.329s
time bash 4_proc_cmd.sh
bash 5_circos_cmd.sh
```

### Atha: partition sequences

```bash
cd ~/data/alignment/self

mkdir -p ~/data/alignment/self/arabidopsis_parted/Genomes
perl ~/Scripts/egaz/part_seq.pl \
    -i ~/data/alignment/Ensembl/Atha \
    -o ~/data/alignment/self/arabidopsis_parted/Genomes/Atha \
    --chunk 10010000 --overlap 10000

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --parted \
    --length 1000  \
    --norm \
    --name arabidopsis_parted \
    --parallel 8 \
    -t Atha

bash arabidopsis_parted/1_real_chr.sh
# real    21m10.875s
time bash arabidopsis_parted/3_self_cmd.sh
# real    9m33.086s
time bash arabidopsis_parted/4_proc_cmd.sh
bash arabidopsis_parted/5_circos_cmd.sh
```

# All plants

## Plants: full chromosomes

```bash
cd ~/data/alignment/self

perl ~/Scripts/egaz/self_batch.pl \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
    --norm \
    --name plants \
    --parallel 12 \
    -q Alyr \
    -q Bdis \
    -q OsatJap \
    -q Sbic \
    -t Atha

cd ~/data/alignment/self/plants

bash 1_real_chr.sh
# real    658m10.559s
time bash 3_self_cmd.sh
# real    971m49.779s
time bash 4_proc_cmd.sh
bash 5_circos_cmd.sh
```

## Plants: partitioned chromosomes

`Bole` contains exact matched pieces with copy number large than one thousand.

It turns out that poor assemblies tend to have this phenomenon.

To speed up processing, use partitioned sequences in `3_self_cmd.sh` and `--discard 50` in
`4_proc_cmd.sh`.

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
    --working_dir ~/data/alignment/self \
    -c ~/data/alignment/self/ensembl_taxon.csv \
    --length 1000  \
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

