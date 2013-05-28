# for f in `find . -name "real_chr.sh" `;do   echo sh $f ;done | sort > taxon_chr_length.sh
# sh taxon_chr_length.sh > merge_taxon_chr_length.sh
# sh merge_taxon_chr_length.sh

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/file-rm.sh ; \
    echo sh $d/pair_cmd.sh ; \
    echo sh $d/multi_cmd.sh ; \
    echo ; \
done  > runall.sh

for d in `find $PWD -mindepth 1 -maxdepth 1 -type d | sort `;do \
    echo sh $d/real_chr.sh ; \
    echo ; \
done  > run_chr.sh
sh run_chr.sh > run_real_chr.sh
sh run_real_chr.sh


# perl d:/wq/Scripts/tool/replace.pl -d d:/wq/Scripts/alignDB/bac -p "cmd.bat" -f /home/wangq -r d:/wq

