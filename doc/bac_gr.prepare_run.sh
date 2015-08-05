#!/bin/bash
cd ~/data/bacteria/process_new

# cat ~/Scripts/withncbi/doc/bac_gr.cmd.txt | grep perl | parallel --jobs 12 --keep-order  2>&1 | tee log_bac_gr.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/prepare.sh ; \
    echo ; \
done  > run_prepare.sh
cat run_prepare.sh | grep . | parallel -j 15 -k 2>&1 | tee log_prepare.txt
    
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/1_real_chr.sh ; \
    echo sh $d/2_file_rm.sh ; \
    echo sh $d/3_pair_cmd.sh ; \
    echo sh $d/4_rawphylo.sh ; \
    echo sh $d/5_multi_cmd.sh ; \
    echo sh $d/6_multi_db_only.sh ; \
    echo ; \
done  > runall.sh

# 1_real_chr
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/1_real_chr.sh ; \
    echo ; \
done  > run_1_real_chr.sh
cat run_1_real_chr.sh | grep . | parallel -j 2 2>&1 | tee log_1_real_chr.txt

# RepeatMasker
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/2_file_rm.sh ; \
    echo ; \
done  > run_2_file_rm.sh
cat run_2_file_rm.sh | grep . | parallel -j 2 2>&1 | tee log_2_file_rm.txt

# pair
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/3_pair_cmd.sh ; \
    echo ; \
done  > run_3_pair_cmd.sh
cat run_3_pair_cmd.sh | grep . | parallel -j 4 2>&1 | tee log_3_pair_cmd.txt

# rawphylo
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/4_rawphylo.sh ; \
    echo ; \
done  > run_4_rawphylo.sh
cat run_4_rawphylo.sh | grep . | parallel -j 4 2>&1 | tee log_4_rawphylo.txt

# multi cmd
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/5_multi_cmd.sh ; \
    echo ; \
done  > run_5_multi_cmd.sh
cat run_5_multi_cmd.sh | grep . | parallel -j 4 2>&1 | tee log_5_multi_cmd.txt

# multi db
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/6_multi_db_only.sh ; \
    echo ; \
done  > run_6_multi_db_only.sh
cat run_6_multi_db_only.sh | grep . | parallel -j 4 2>&1 | tee log_6_multi_db_only.txt

# clean
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_raw" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "rawphylo" | xargs rm -fr 

find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy" | xargs rm 
find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy.reduced" | xargs rm 

# split -d --lines=100 runall.sh runall.
# sh runall.00

# split -d --lines=100 run_multi_db_only.sh run_multi_db_only.

