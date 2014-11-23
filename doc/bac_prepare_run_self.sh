#!/bin/bash
cd ~/data/bacteria/process_self

# cat ~/Scripts/alignDB/bac/species_gr_self_cmd.txt | grep perl | parallel --jobs 12 --keep-order  2>&1 | tee log_gr.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/prepare.sh ; \
    echo ; \
done  > run_prepare.sh
cat run_prepare.sh | grep . | parallel -j 15 -k 2>&1 | tee log_prepare.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/real_chr.sh ; \
    echo ; \
done  > run_chr.sh
sh run_chr.sh > run_real_chr.sh
sh run_real_chr.sh

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/file-rm.sh ; \
    echo sh $d/self_cmd.sh ; \
    echo sh $d/proc_cmd.sh ; \
    echo sh $d/circos_cmd.sh ; \
    echo ; \
done  > runall.sh

# RepeatMasker
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/file-rm.sh ; \
    echo ; \
done  > run_file_rm.sh
cat run_file_rm.sh | grep . | parallel -j 2 2>&1 | tee log_file_rm.txt

# generate new circos.conf
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/redo_prepare.sh ; \
    echo ; \
done  > run_redo_prepare.sh
cat run_redo_prepare.sh | grep . | parallel -j 15

# lastz
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/self_cmd.sh ; \
    echo ; \
done  > run_self_cmd.sh
cat run_self_cmd.sh | grep . | parallel -j 15 2>&1 | tee log_self_cmd.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/proc_cmd.sh ; \
    echo ; \
done  > run_proc_cmd.sh
cat run_proc_cmd.sh | grep . | parallel -j 15 2>&1 | tee log_proc.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/circos_cmd.sh ; \
    echo ; \
done  > run_circos_cmd.sh
cat run_circos_cmd.sh | grep . | parallel -j 15 2>&1 | tee log_circos.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/feature_cmd.sh ; \
    echo ; \
done  > run_feature_cmd.sh
cat run_feature_cmd.sh | grep . | parallel -j 15 2>&1 | tee log_feature.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/pair_stat.sh ; \
    echo ; \
done  > run_pair_stat.sh
cat run_pair_stat.sh | grep . | parallel -j 4 2>&1 | tee log_pair_stat.txt

# clean
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_raw" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "rawphylo" | xargs rm -fr 

find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy" | xargs rm 
find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy.reduced" | xargs rm 

# split -d --lines=100 runall.sh runall.
# sh runall.00

# split -d --lines=100 run_multi_db_only.sh run_multi_db_only.

### return to state before run_proc_cmd.sh
find $PWD -type f -name *xlsx | xargs rm
find $PWD -type f -name *yml | xargs rm
find $PWD -type d -name *_proc | xargs rm -fr
find $PWD -type d -name *_result | xargs rm -fr
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/redo_prepare.sh ; \
    echo ; \
done  > run_redo_prepare.sh
cat run_redo_prepare.sh | grep . | parallel -j 15

### summary
mkdir ~/data/bacteria/process_self_summary
cd ~/data/bacteria/process_self_summary

# $ARGV is the name of the current file
# $ARGV[0] *isn't* !
find  ~/data/bacteria/process_self -type f -name "*.1000.csv" | sort \
    | xargs perl -nle'BEGIN{print q{species,strain_id,group,chr,length,size,coverage}}; /^u\,all/ and $ARGV =~ /\/(\w+)\/(\d+)_/ and print qq{$1,$2,$_}' \
    > bacteria.paralog.union.csv

find  ~/data/bacteria/process_self -type f -name "*.cc.csv" | sort \
    | xargs perl -nle'BEGIN{print q{species,strain_id,copy,chr,length,size,coverage,percent,count,piece,avg_size}}; /\,sum/ and $ARGV =~ /\/(\w+)\/(\d+)_/ and print qq{$1,$2,$_}' \
    > bacteria.paralog.copies.csv

find  ~/data/bacteria/process_self -type f -name "*.feature.csv" | sort \
    | xargs perl -nle'BEGIN{print q{species,strain_id,feature,chr,length,size,coverage}}; /\,all/ and $ARGV =~ /\/(\w+)\/(\d+)_/ and print qq{$1,$2,$_}' \
    > bacteria.feature.csv

find  ~/data/bacteria/process_self -type f -name "*.feature.copies.csv" | sort \
    | xargs perl -nle'BEGIN{print q{species,strain_id,feature,copy,chr,length,size,coverage}}; /\,all/ and $ARGV =~ /\/(\w+)\/(\d+)_/ and print qq{$1,$2,$_}' \
    > bacteria.feature.copies.csv

find  ~/data/bacteria/process_self -type f -name "info.csv" | sort \
    | xargs perl -nle'BEGIN{print q{strain,strain_id,species,species_id,genus,genus_id,family,family_id,order,order_id}}; /^strain\,/ and next; print $_' \
    > bacteria.species.info.csv

find  ~/data/bacteria/process_self -type f -name "chr.sizes" | sort \
    | xargs perl -nle'BEGIN{print q{species,strain_id,accession,length}}; $_ =~ s/\t/\,/; $ARGV =~ /process_self\/(\w+)\/(\d+)\//; print qq{$1,$2,$_}' \
    > bacteria.species.genome.csv

find  ~/data/bacteria/process_self -type f -name "chr.sizes" | sort \
    | perl -nle'/process_self/ or next; @f = grep {defined } split /\//; $f[-1] =~ s/\.fa//; print qq{$f[-3],$f[-2],$f[-1]}'
    
     \
    | sort | uniq | xargs mkdir 



# count feature
find  ~/data/bacteria/process_self -type f -name "*.gff" \
    | xargs perl -anle'/^#/ and next; print $F[2]' \
    | sort | uniq -c \
    > bacteria.feature.count

find  ~/data/bacteria/process_self -type f -name "*.gff" \
    | xargs perl -anle'/^#/ and next; $F[2] eq 'dispersed_repeat' and print qq{$F[2]\t$F[8]}' \
    | sort | uniq -c

find  ~/data/bacteria/process_self -type f -name "*.gff" \
    | xargs perl -anle'/^#/ and next; $F[2] eq 'region' and $F[8] =~ /note\=(\w+)/i and print qq{$F[2]\t$1}' \
    | sort | uniq -c \
    | perl -anle'$F[0] > 10 and print' \
    > bacteria.feature.region.note.count

find  ~/data/bacteria/process_self -type f -name "*.gff" \
    | xargs perl -anle'/^#/ and next; $F[2] eq 'region' and $F[8] =~ /gbkey\=(\w+)/i and print qq{$F[2]\t$1}' \
    | sort | uniq -c \
    > bacteria.feature.region.gbkey.count

find  ~/data/bacteria/process_self -type f -name "*.gff" \
    | xargs perl -anle'/^#/ and next; $F[2] eq 'repeat_region' and print qq{$F[2]\t$F[8]}' \
    | sort | uniq -c \
    > bacteria.feature.repeat_region.count

find  ~/data/bacteria/process_self -type f -name "*.gff" \
    | xargs perl -anle'/^#/ and next; $F[2] eq 'dispersed_repeat' and $F[8] !~ /rna/i and print qq{$F[2]\t$F[8]}' \
    | sort | uniq -c \
    > bacteria.feature.dispersed_repeat.count


### xlsx
# mkdir ~/data/bacteria/process_self_xlsx
cd ~/data/bacteria/process_self_xlsx

find  ~/data/bacteria/process_self -type f -name "*basicstat.xlsx" | sort \
    | xargs -i cp {} .

find  ~/data/bacteria/process_self -type f -name "*_paralog.common.xlsx" | sort \
    | xargs -i cp {} .

find  ~/data/bacteria/process_self -type f -name "*_paralog.gc.xlsx" | sort \
    | xargs -i cp {} .

find  ~/data/bacteria/process_self -type f -name "*_paralog.mvar.xlsx" | sort \
    | xargs -i cp {} .


### circos
# mkdir ~/data/bacteria/process_self_circos
cd ~/data/bacteria/process_self_circos

find  ~/data/bacteria/process_self -type f -name "*.png" | sort \
    | perl -nle'/process_self\/(\w+)\/.+\/(\d+\.circos\.png)/ and print qq{$1}' \
    | sort | uniq | xargs mkdir 

find  ~/data/bacteria/process_self -type f -name "*.png" | sort \
    | perl -nle'/process_self\/(\w+)\/.+\/(\d+\.circos\.png)/ and system qq{cp $_ $1}'


### species
# mkdir ~/data/bacteria/process_self_tar
cd ~/data/bacteria/process_self_tar

find  ~/data/bacteria/process_self -maxdepth 1 -mindepth 1 -type d | xargs -i basename {} | sort \
    | parallel -j 6 tar zcvf {}.tar.gz ~/data/bacteria/process_self/{}

### phylo
# mkdir ~/data/bacteria/process_self_phylo
cd ~/data/bacteria/process_self_phylo

find  ~/data/bacteria/process_gr -type f -name "RAxML_bipartitions.*" | sort \
    | perl -nle'/RAxML_bipartitions\.(\w+)/ and print qq{$1}' \
    | sort | uniq | xargs mkdir 

find  ~/data/bacteria/process_gr -type f -name "RAxML_bipartitions.*" | sort \
    | perl -nle'/RAxML_bipartitions\.(\w+)/ and system qq{cp $_ $1\.nwk}'


# gr xlsx
# mkdir ~/data/bacteria/process_gr_xlsx
cd ~/data/bacteria/process_gr_xlsx
find  ~/data/bacteria/process_gr -type f -name "*.common.chart.xlsx" | sort \
    | xargs -i cp {} .

find  ~/data/bacteria/process_gr -type f -name "*.multi.xlsx" | sort \
    | xargs -i cp {} .

find  ~/data/bacteria/process_gr -type f -name "*.gc.xlsx" | sort \
    | xargs -i cp {} .

find  ~/data/bacteria/process_gr -type f -name "*.mvar.xlsx" | sort \
    | xargs -i cp {} .

# table.txt
# mkdir ~/data/bacteria/process_gr_table_txt
cd ~/data/bacteria/process_gr_table_txt
find  ~/data/bacteria/process_gr -type f -name "table.txt" | sort \
    | perl -nle'/process_gr\/(\w+)/ and system qq{cp $_ $1\.table\.txt}'

