# Self-aligning steps for each groups


[TOC levels=1-3]: # ""

- [Self-aligning steps for each groups](#self-aligning-steps-for-each-groups)
- [Taxonomy for Ensembl species](#taxonomy-for-ensembl-species)
  - [Plants](#plants)
- [Model organisms](#model-organisms)
  - [E. coli K-12 MG1655](#e-coli-k-12-mg1655)
  - [Yeast S288c](#yeast-s288c)
  - [Dmel](#dmel)
  - [Cele](#cele)
  - [Ddis](#ddis)
  - [Human](#human)
  - [Mouse](#mouse)
- [All plants](#all-plants)
  - [Arabidopsis](#arabidopsis)
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

## E. coli K-12 MG1655

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
spanr gff --tag CDS chr.gff -o cds.yml
spanr gff U00096.rm.gff -o repeat.yml
spanr merge repeat.yml cds.yml -o anno.yml

rm repeat.yml cds.yml U00096.rm.gff U00096.rm.out
rm -fr MG1655

faops masked *.fa |
    spanr cover stdin |
    spanr stat --all chr.sizes stdin
#chrLength,size,coverage
#4641652,21767,0.0047

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

faops masked ~/data/alignment/Ensembl/S288c/*.fa |
    spanr cover stdin |
    spanr stat --all ~/data/alignment/Ensembl/S288c/chr.sizes stdin
#chrLength,size,coverage
#12071326,729220,0.0604

egaz template \
    ~/data/alignment/Ensembl/S288c/ \
    --self -o yeast \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 8 -v

bash yeast/1_self.sh
bash yeast/3_proc.sh
bash yeast/4_circos.sh

cat yeast/Results/S288c/S288c.cover.csv |
    grep "^all"
#all,12071326,695510,0.0576

```

## Dmel

```bash
cd ~/data/alignment/self

faops masked ~/data/alignment/Ensembl/Dmel/*.fa |
    spanr cover stdin |
    spanr stat --all ~/data/alignment/Ensembl/Dmel/chr.sizes stdin
#chrLength,size,coverage
#132532477,8124132,0.0613

egaz template \
    ~/data/alignment/Ensembl/Dmel/ \
    --self -o fly \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 16 -v

bash fly/1_self.sh
bash fly/3_proc.sh
bash fly/4_circos.sh

cat fly/Results/Dmel/Dmel.cover.csv |
    grep "^all"
#all,132532477,13971696,0.1054

```

## Cele

```bash
cd ~/data/alignment/self

faops masked ~/data/alignment/Ensembl/Cele/*.fa |
    spanr cover stdin |
    spanr stat --all ~/data/alignment/Ensembl/Cele/chr.sizes stdin
#chrLength,size,coverage
#100272607,22008206,0.2195

egaz template \
    ~/data/alignment/Ensembl/Cele/ \
    --self -o worm \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 16 -v

bash worm/1_self.sh
bash worm/3_proc.sh
bash worm/4_circos.sh

cat worm/Results/Cele/Cele.cover.csv |
    grep "^all"
#all,100272607,4702767,0.0469

```

## Ddis

* .dna_sm.toplevel.fa.gz is unmasked

```bash
rm -fr  ~/data/alignment/Ensembl/Ddis/

egaz prepseq \
    --repeatmasker '--gff --parallel 8' --min 50000 -v \
    ~/data/ensembl98/fasta/dictyostelium_discoideum/dna/Dictyostelium_discoideum.dicty_2.7.dna_sm.toplevel.fa.gz \
    -o ~/data/alignment/Ensembl/Ddis/

cd ~/data/alignment/Ensembl/Ddis/

find ~/data/ensembl98/gff3/dictyostelium_discoideum/ -name "*.gff3.gz" |
    grep -v "abinitio.gff3" |
    grep -v "chr.gff3" |
    xargs gzip -d -c > chr.gff
    
# create anno.yml
spanr gff --tag CDS chr.gff -o cds.yml
spanr gff *.rm.gff -o repeat.yml
spanr merge repeat.yml cds.yml -o anno.yml

rm -f repeat.yml cds.yml

```

```bash
cd ~/data/alignment/self

faops masked ~/data/alignment/Ensembl/Ddis/*.fa |
    spanr cover stdin |
    spanr stat --all ~/data/alignment/Ensembl/Ddis/chr.sizes stdin
#chrLength,size,coverage
#33943072,7442427,0.2193

egaz template \
    ~/data/alignment/Ensembl/Ddis/ \
    --self -o dicty \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 16 -v

time bash dicty/1_self.sh
#real    5m31.495s
#user    68m13.782s
#sys     1m22.881s

bash dicty/3_proc.sh
bash dicty/4_circos.sh

cat dicty/Results/Ddis/Ddis.cover.csv |
    grep "^all"
#all,33943072,4568080,0.1346

```

## Human

```bash
cd ~/data/alignment/self

egaz template \
    ~/data/alignment/Ensembl/Human/ \
    --self -o human \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 16 -v

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
    --circos --parallel 16 -v

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

Comparison: **9.58% vs 8.65%**. Use full chromosomes if the running time is acceptable.

```bash
cd ~/data/alignment/self

faops masked ~/data/alignment/Ensembl/Atha/*.fa |
    spanr cover stdin |
    spanr stat --all ~/data/alignment/Ensembl/Atha/chr.sizes stdin
#chrLength,size,coverage
#119146348,28109149,0.2359

# full chromosomes
egaz template \
    ~/data/alignment/Ensembl/Atha/ \
    --self -o arabidopsis \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 16 -v

time bash arabidopsis/1_self.sh
#real    14m46.549s
#user    166m24.298s
#sys     1m38.050s

time bash arabidopsis/3_proc.sh
#real    7m22.867s
#user    27m4.939s
#sys     14m4.978s

bash arabidopsis/4_circos.sh

cat arabidopsis/Results/Atha/Atha.cover.csv |
    grep "^all"
#all,119146348,11413522,0.0958

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
    --circos --partition --parallel 16 -v

time bash arabidopsis_par/1_self.sh
#real    11m27.779s
#user    168m57.061s
#sys     4m52.982s

time bash arabidopsis_par/3_proc.sh
#real    6m11.381s
#user    23m23.330s
#sys     12m57.028s

bash arabidopsis_par/4_circos.sh

cat arabidopsis_par/Results/Atha/Atha.cover.csv |
    grep "^all"
#all,119146348,10305900,0.0865

```

## Plants: full chromosomes

```bash
cd ~/data/alignment/self

for name in Atha Alyr OsatJap Sbic; do
    echo ${name}
    faops masked ~/data/alignment/Ensembl/${name}/*.fa |
        spanr cover stdin |
        spanr stat --all ~/data/alignment/Ensembl/${name}/chr.sizes stdin
done
#Atha
#chrLength,size,coverage
#119146348,28109149,0.2359
#Alyr
#chrLength,size,coverage
#194182311,82250191,0.4236
#OsatJap
#chrLength,size,coverage
#373245519,190136060,0.5094
#Sbic
#chrLength,size,coverage
#683645045,468096621,0.6847

egaz template \
    ~/data/alignment/Ensembl/Atha/ \
    ~/data/alignment/Ensembl/Alyr/ \
    ~/data/alignment/Ensembl/OsatJap/ \
    ~/data/alignment/Ensembl/Sbic/ \
    --self -o plants \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --parallel 16 -v

time bash plants/1_self.sh
#real    133m47.288s
#user    1387m27.563s
#sys     16m55.046s

time bash plants/3_proc.sh

bash plants/4_circos.sh

```

## Plants: partitioned chromosomes

```bash
cd ~/data/alignment/self

for name in Atha Alyr OsatJap Sbic Mtru Gmax Bole Brap Vvin Slyc Stub Macu Sita Bdis; do
    echo ${name}
    faops masked ~/data/alignment/Ensembl/${name}/*.fa |
        spanr cover stdin |
        spanr stat --all ~/data/alignment/Ensembl/${name}/chr.sizes stdin
done
#Atha
#chrLength,size,coverage
#119146348,28109149,0.2359
#Alyr
#chrLength,size,coverage
#194182311,82250191,0.4236
#OsatJap
#chrLength,size,coverage
#373245519,190136060,0.5094
#Sbic
#chrLength,size,coverage
#683645045,468096621,0.6847
#Mtru
#chrLength,size,coverage
#384466993,145484881,0.3784
#Gmax
#chrLength,size,coverage
#949183385,478117318,0.5037
#Bole
#chrLength,size,coverage
#446885882,133157067,0.2980
#Brap
#chrLength,size,coverage
#256240462,55286796,0.2158
#Vvin
#chrLength,size,coverage
#426176009,184418981,0.4327
#Slyc
#chrLength,size,coverage
#807224664,493504368,0.6114
#Stub
#chrLength,size,coverage
#810654046,527639148,0.6509
#Macu
#chrLength,size,coverage
#331812599,88618121,0.2671
#Sita
#chrLength,size,coverage
#401296418,107338036,0.2675
#Bdis
#chrLength,size,coverage
#271067295,11899938,0.0439

for name in Atha Alyr OsatJap Sbic Mtru Gmax Bole Brap Vvin Slyc Stub Macu Sita Bdis; do
    echo "==> ${name}"
    find ~/data/alignment/Ensembl/${name}/ -type f -name "*.fa" |
        sort |
        parallel --no-run-if-empty --linebuffer -k -j 8 '
            echo >&2 {}
            egaz partition {} --chunk 10010000 --overlap 10000
        '
done

egaz template \
    ~/data/alignment/Ensembl/Atha/ \
    ~/data/alignment/Ensembl/Alyr/ \
    ~/data/alignment/Ensembl/OsatJap/ \
    ~/data/alignment/Ensembl/Sbic/ \
    ~/data/alignment/Ensembl/Mtru/ \
    ~/data/alignment/Ensembl/Gmax/ \
    ~/data/alignment/Ensembl/Bole/ \
    ~/data/alignment/Ensembl/Brap/ \
    ~/data/alignment/Ensembl/Vvin/ \
    ~/data/alignment/Ensembl/Slyc/ \
    ~/data/alignment/Ensembl/Stub/ \
    ~/data/alignment/Ensembl/Bdis/ \
    --self -o plants_par \
    --taxon ~/data/alignment/self/ensembl_taxon.csv \
    --circos --partition --parallel 16 -v

time bash plants_par/1_self.sh
#real    3780m34.746s
#user    39108m41.641s
#sys     1656m31.650s

time bash plants_par/3_proc.sh
#real    1177m51.175s
#user    2609m6.773s
#sys     1307m16.303s

bash plants_par/4_circos.sh

```

* Results of `egaz template`

```text
$ du -hs ~/data/alignment/self/plants_par/Pairwise/*
332M    /home/wangq/data/alignment/self/plants_par/Pairwise/AlyrvsSelf
111M    /home/wangq/data/alignment/self/plants_par/Pairwise/AthavsSelf
25G     /home/wangq/data/alignment/self/plants_par/Pairwise/BdisvsSelf
25G     /home/wangq/data/alignment/self/plants_par/Pairwise/BolevsSelf
3.6G    /home/wangq/data/alignment/self/plants_par/Pairwise/BrapvsSelf
8.6G    /home/wangq/data/alignment/self/plants_par/Pairwise/GmaxvsSelf
2.7G    /home/wangq/data/alignment/self/plants_par/Pairwise/MtruvsSelf
263M    /home/wangq/data/alignment/self/plants_par/Pairwise/OsatJapvsSelf
1.6G    /home/wangq/data/alignment/self/plants_par/Pairwise/SbicvsSelf
8.5G    /home/wangq/data/alignment/self/plants_par/Pairwise/SlycvsSelf
4.0G    /home/wangq/data/alignment/self/plants_par/Pairwise/StubvsSelf
14G     /home/wangq/data/alignment/self/plants_par/Pairwise/VvinvsSelf

```

* Results of the old `egaz/self_batch.pl`

```text
$ du -hs ~/data/alignment/self/plants_parted/Pairwise/*
890M    /home/wangq/data/alignment/self/plants_parted/Pairwise/Alyrvsselfalign
94M     /home/wangq/data/alignment/self/plants_parted/Pairwise/Athavsselfalign
6.2G    /home/wangq/data/alignment/self/plants_parted/Pairwise/Bdisvsselfalign
17G     /home/wangq/data/alignment/self/plants_parted/Pairwise/Bolevsselfalign
2.2G    /home/wangq/data/alignment/self/plants_parted/Pairwise/Brapvsselfalign
42G     /home/wangq/data/alignment/self/plants_parted/Pairwise/Gmaxvsselfalign
4.4G    /home/wangq/data/alignment/self/plants_parted/Pairwise/Mtruvsselfalign
845M    /home/wangq/data/alignment/self/plants_parted/Pairwise/OsatJapvsselfalign
103G    /home/wangq/data/alignment/self/plants_parted/Pairwise/Slycvsselfalign
15G     /home/wangq/data/alignment/self/plants_parted/Pairwise/Stubvsselfalign
9.3G    /home/wangq/data/alignment/self/plants_parted/Pairwise/Vvinvsselfalign

```

