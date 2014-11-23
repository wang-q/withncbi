#!/bin/bash
cd ~/data/bacteria/process_withncbi

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
    echo sh $d/pair_cmd.sh ; \
    echo sh $d/rawphylo.sh ; \
    echo sh $d/multi_cmd.sh ; \
    echo ; \
done  > runall.sh

# RepeatMasker
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/file-rm.sh ; \
    echo ; \
done  > run_file_rm.sh
cat run_file_rm.sh | grep . | parallel -j 2 2>&1 | tee log_file_rm.txt
# pair
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/pair_cmd.sh ; \
    echo ; \
done  > run_pair_cmd.sh
cat run_pair_cmd.sh | grep . | parallel -j 15 2>&1 | tee log_pair_cmd.txt
for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/rawphylo.sh ; \
    echo ; \
done  > run_rawphylo.sh
cat run_rawphylo.sh | grep . | parallel -j 15 2>&1 | tee log_rawphylo.txt

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/multi_cmd.sh ; \
    echo ; \
done  > run_multi_cmd.sh
cat run_multi_cmd.sh | grep . | parallel -j 4 2>&1 | tee log_multi_cmd.txt
# clean
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_raw" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "rawphylo" | xargs rm -fr 

find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy" | xargs rm 
find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy.reduced" | xargs rm 

# split -d --lines=100 runall.sh runall.
# sh runall.00

# split -d --lines=100 run_multi_db_only.sh run_multi_db_only.

