
# sh ~/Scripts/alignDB/bac/species_gr_woo_cmd.txt 2>&1 | tee log_gr.txt

# cd ~/data/bacteria/process_gr



for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/prepare.sh ; \
    echo ; \
done  > run_prepare.sh

# sh run_prepare.sh 2>&1 | tee log_prepare.txt

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

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/multi_db_only.sh ; \
    echo ; \
done  > run_multi_db_only.sh

# clean
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_raw" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "*_fasta" | xargs rm -fr 
find $PWD -mindepth 1 -maxdepth 2 -type d -name "rawphylo" | xargs rm -fr 

find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy" | xargs rm 
find $PWD -mindepth 1 -maxdepth 3 -type f -name "*.phy.reduced" | xargs rm 

# split --lines=100 runall.sh runall.
# sh runall.aa

# split -d --lines=100 run_multi_db_only.sh run_multi_db_only.

