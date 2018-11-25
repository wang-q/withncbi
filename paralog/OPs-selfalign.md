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
    --id    9606   --name 9606=Human    \
    --id    10090  --name 10090=Mouse   \
    --id    3702   --name 3702=Atha     \
    --id    59689  --name 59689=Alyr    \
    --id    3880   --name 3880=Mtru     \
    --id    3847   --name 3847=Gmax     \
    --id    29760  --name 29760=Vvin    \
    --id    39947  --name 39947=OsatJap \
    --id    3712   --name 3712=Bole     \
    --id    3711   --name 3711=Brap     \
    --id    4081   --name 4081=Slyc     \
    --id    4113   --name 4113=Stub     \
    --id    4555   --name 4555=Sita     \
    --id    4558   --name 4558=Sbic     \
    --id    15368  --name 15368=Bdis    \
    --id    4641   --name 4641=Macu     \
    --id    7227   --name 7227=Dmel     \
    --id    6239   --name 6239=Cele     \
    --id    352472 --name 352472=Ddis   \
    --file ensembl_taxon.csv            \
    --entrez

```

## Plants

| Name                       | Classification | Taxon ID | Used |
|:---------------------------|:---------------|:---------|:----:|
| Aegilops tauschii          | Liliopsida     | 37682    |      |
| Amborella trichopoda       | Amborellales   | 13333    |  o   |
| Arabidopsis lyrata         | eudicotyledons | 81972    |  o   |
| Arabidopsis thaliana       | eudicotyledons | 3702     |  o   |
| Brachypodium distachyon    | Liliopsida     | 15368    |  o   |
| Brassica oleracea          | eudicotyledons | 109376   |      |
| Brassica rapa              | eudicotyledons | 51351    |  o   |
| Chlamydomonas reinhardtii  | Chlorophyta    | 3055     |      |
| Cyanidioschyzon merolae    | Rhodophyta     | 280699   |      |
| Glycine max                | eudicotyledons | 3847     |  o   |
| Hordeum vulgare            | Liliopsida     | 112509   |      |
| Leersia perrieri           | Liliopsida     | 77586    |      |
| Medicago truncatula        | eudicotyledons | 3880     |  o   |
| Musa acuminata             | Liliopsida     | 214687   |  o   |
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
| Oryza sativa Japonica      | Liliopsida     | 39947    |  o   |
| Ostreococcus lucimarinus   | Chlorophyta    | 436017   |      |
| Physcomitrella patens      | Bryophyta      | 3218     |      |
| Populus trichocarpa        | eudicotyledons | 3694     |      |
| Prunus persica             | eudicotyledons | 3760     |      |
| Selaginella moellendorffii | Lycopodiophyta | 88036    |      |
| Setaria italica            | Liliopsida     | 4555     |  o   |
| Solanum lycopersicum       | eudicotyledons | 4081     |  o   |
| Solanum tuberosum          | eudicotyledons | 4113     |  o   |
| Sorghum bicolor            | Liliopsida     | 4558     |  o   |
| Theobroma cacao            | eudicotyledons | 3641     |      |
| Triticum aestivum          | Liliopsida     | 4565     |      |
| Triticum urartu            | Liliopsida     | 4572     |      |
| Vitis vinifera             | eudicotyledons | 29760    |  o   |
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

perl ~/Scripts/withncbi/taxon/batch_get_seq.pl -f name_seq.csv

egaz prepseq MG1655/U00096.fa -o . --repeatmasker '--gff --parallel 8' -v

mv MG1655/U00096.gff .
mv U00096.gff chr.gff

# create anno.yml
runlist gff --tag CDS --remove chr.gff -o cds.yml
runlist gff --remove U00096.rm.gff -o repeat.yml
runlist merge repeat.yml cds.yml -o anno.yml

rm repeat.yml cds.yml U00096.rm.gff U00096.rm.out
rm -fr MG1655

```

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/MG1655 \
    --self -o ecoli \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --aligndb --parallel 8 -v

bash ecoli/1_self.sh
bash ecoli/3_proc.sh
bash ecoli/4_circos.sh
bash ecoli/6_chr_length.sh
bash ecoli/7_self_aligndb.sh
bash ecoli/9_pack_up.sh

```

## Yeast S288c

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/S288c/ \
    --self -o yeast \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

bash yeast/1_self.sh
bash yeast/3_proc.sh
bash yeast/4_circos.sh

```

## Dmel

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/Dmel/ \
    --self -o fly \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

time bash fly/1_self.sh
#real    20m4.269s
#user    139m8.359s
#sys     1m42.665s

time bash fly/3_proc.sh
#real    6m16.450s
#user    16m32.193s
#sys     5m26.875s

bash fly/4_circos.sh

```

## Cele

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/Cele/ \
    --self -o worm \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

bash worm/1_self.sh
bash worm/3_proc.sh
bash worm/4_circos.sh

```

## Ddis

* .dna_sm.toplevel.fa.gz is unmasked

```bash
rm -fr  ~/data/alignment/Ensembl/Ddis/

egaz prepseq \
    --repeatmasker '--gff --parallel 8' --min 50000 -v \
    ~/data/ensembl94/fasta/dictyostelium_discoideum/dna/Dictyostelium_discoideum.dicty_2.7.dna_sm.toplevel.fa.gz \
    -o ~/data/alignment/Ensembl/Ddis/

cd ~/data/alignment/Ensembl/Ddis/

find ~/data/ensembl94/gff3/dictyostelium_discoideum/ -name "*.gff3.gz" |
    grep -v "abinitio.gff3" |
    grep -v "chr.gff3" |
    xargs gzip -d -c > chr.gff
runlist gff --tag CDS --remove chr.gff -o cds.yml

runlist gff --remove \
    *.rm.gff \
    -o repeat.yml

runlist merge \
    cds.yml repeat.yml \
    -o anno.yml

rm -f repeat.yml cds.yml

```

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/Ddis/ \
    --self -o dicty \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

time bash dicty/1_self.sh
#real    7m32.264s
#user    54m42.929s
#sys     0m30.717s

bash dicty/3_proc.sh
bash dicty/4_circos.sh

```

## Human

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/Human/ \
    --self -o human \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

time bash human/1_self.sh
#real    3158m12.839s
#user    24898m37.946s
#sys     138m1.231s

time bash human/3_proc.sh
#real    458m46.491s
#user    1330m6.094s
#sys     93m22.268s

bash human/4_circos.sh

```

## Mouse

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/Mouse/ \
    --self -o mouse \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

time bash mouse/1_self.sh
#real    2948m46.920s
#user    23114m31.078s
#sys     108m47.712s

time bash mouse/3_proc.sh
#real    1012m29.622s
#user    2487m21.989s
#sys     193m31.933s

bash mouse/4_circos.sh

```

# All plants

## Arabidopsis

Comparison: **9.69% vs 9.34%**. Use full chromosomes if the running time is acceptable.

```bash
cd ~/data/alignment/self

# full chromosomes
egaz template \
    ~/data/alignment/Ensembl/Atha/ \
    --self -o arabidopsis \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

time bash arabidopsis/1_self.sh
#real    17m35.952s
#user    121m5.877s
#sys     0m23.514s

time bash arabidopsis/3_proc.sh
#real    8m59.367s
#user    23m0.583s
#sys     14m46.156s

bash arabidopsis/4_circos.sh

# partitioned chromosomes
find ~/data/alignment/Ensembl/Atha/ -type f -name "*.fa" |
    sort |
    parallel --no-run-if-empty --linebuffer -k -j 8 '
        echo >&2 {}
        egaz partition {} --chunk 10010000 --overlap 10000
    '

egaz template \
    ~/data/alignment/Ensembl/Atha/ \
    --self -o arabidopsis_par \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --partition --parallel 8 -v

time bash arabidopsis_par/1_self.sh
#real    16m47.681s
#user    128m26.793s
#sys     0m56.079s

time bash arabidopsis_par/3_proc.sh
#real    9m5.014s
#user    21m50.142s
#sys     13m54.494s

bash arabidopsis_par/4_circos.sh

```

## Plants: full chromosomes

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/Atha/ \
    ~/data/alignment/Ensembl/Alyr/ \
    ~/data/alignment/Ensembl/OsatJap/ \
    ~/data/alignment/Ensembl/Sbic/ \
    --self -o plants \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 12 -v

time bash plants/1_self.sh
#real    133m47.288s
#user    1387m27.563s
#sys     16m55.046s

time bash plants/3_proc.sh
bash plants/4_circos.sh

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

