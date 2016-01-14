# Self-aligning steps for each groups

## Taxonomy for Ensembl species

```bash
mkdir -p ~/data/alignment/self
cd ~/data/alignment/self

perl       ~/Scripts/withncbi/taxon/strain_info.pl \
    --file ensembl_taxon.csv                       \
    --id   559292                                  \
    --name 559292=S288c                            \
    --id   3702                                    \
    --name 3702=Atha                               \
    --id   39947                                   \
    --name 39947=OsatJap                           \
    --id   7227                                    \
    --name 7227=Dmel                               \
    --id   9606                                    \
    --name 9606=Human                              \
    --id   6239                                    \
    --name 6239=Cele                               \
    --id   352472                                  \
    --name 352472=Ddis                             \
    --id   10090                                   \
    --name 10090=Mouse                             \
    --entrez

```

## Yeast S288c

All other species got sequences in `OPs-download.md`.

```bash
mkdir -p ~/data/alignment/Ensembl/S288c
cd ~/data/alignment/Ensembl/S288c

find ~/data/ensembl82/fasta/saccharomyces_cerevisiae/dna/ -name "*dna_sm.toplevel*" | xargs gzip -d -c > toplevel.fa
faops count toplevel.fa | perl -aln -e 'next if $F[0] eq 'total'; print $F[0] if $F[1] > 50000; print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05' | uniq > listFile
faops some toplevel.fa listFile toplevel.filtered.fa
faops split-name toplevel.filtered.fa .
rm toplevel.fa toplevel.filtered.fa listFile

mv Mito.fa Mito.fa.skip

# rsync -avP wangq@wq.nju.edu.cn:data/alignment/Ensembl/ ~/data/alignment/Ensembl
```

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
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

sh 1_real_chr.sh
sh 3_self_cmd.sh
sh 4_proc_cmd.sh
sh 5_circos_cmd.sh
```

## Arabidopsis

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
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

sh 1_real_chr.sh
sh 3_self_cmd.sh
time sh 4_proc_cmd.sh
# real    25m5.028s
# user    63m3.210s
# sys     14m29.024s
sh 5_circos_cmd.sh
```

## Rice

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name rice \
    --parallel 8 \
    -t OsatJap

cd ~/data/alignment/self/rice

sh 1_real_chr.sh
sh 3_self_cmd.sh
time sh 4_proc_cmd.sh
# real    280m25.802s
# user    1005m20.730s
# sys     69m0.916s
sh 5_circos_cmd.sh
```

## Fly

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
    --file ~/data/alignment/self/ensembl_taxon.csv \
    --working_dir ~/data/alignment/self \
    --seq_dir ~/data/alignment/Ensembl \
    --length 1000  \
    --use_name \
    --norm \
    --name fly \
    --parallel 8 \
    -t Dmel

cd ~/data/alignment/self/fly

sh 1_real_chr.sh
sh 3_self_cmd.sh
time sh 4_proc_cmd.sh
# real    11m22.152s
# user    31m12.044s
# sys     5m52.670s
sh 5_circos_cmd.sh
```

## Cele

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
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

sh 1_real_chr.sh
time sh 3_self_cmd.sh
# real    26m9.090s
# user    150m18.730s
# sys     0m49.323s
time sh 4_proc_cmd.sh
# real    6m49.232s
# user    12m26.789s
# sys     4m51.724s
sh 5_circos_cmd.sh
```

## Ddis

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
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

sh 1_real_chr.sh
sh 3_self_cmd.sh
time sh 4_proc_cmd.sh
sh 5_circos_cmd.sh
```

## Human

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
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

sh 1_real_chr.sh
time sh 3_self_cmd.sh
# real    4156m4.259s
# user    32775m39.086s
# sys     227m18.093s
time sh 4_proc_cmd.sh
# real    1888m49.720s
# user    7892m55.508s
# sys     220m35.101s
sh 5_circos_cmd.sh
```

## Mouse

```bash
cd ~/data/alignment/self

perl ~/Scripts/withncbi/taxon/strain_bz_self.pl \
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

sh 1_real_chr.sh
time sh 3_self_cmd.sh
time sh 4_proc_cmd.sh
sh 5_circos_cmd.sh
```
