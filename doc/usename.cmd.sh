#----------------------------------------------------------#
# yeast
#----------------------------------------------------------#

# prepare seqs and yeast_ncbi.csv in ~/data/alignment/yeast_genome
S288c       559292 ensembl65
RM11_1a     285006 "NCBI assemble" "10x sanger"
YJM789      307796 "NCBI WGS"      "10x sanger"
Sigma1287b  658763 "NCBI assemble" "45x sanger/solexa"
JAY291      574961 "NCBI WGS"      "12x 454; 58x solexa se; 95x solexa pe"
EC1118      643680 "NCBI WGS"      "17.6x 454; 4.1+1.9x sanger"
Kyokai_no_7 721032 "NCBI WGS"      "9.1x sanger"
T7          929585 "NCBI WGS"      "25.4x 454/sanger"
Spar        226125 "NCBI WGS"      "7x sanger"

# rsync --progress -av wangq@139.162.23.84:/home/wangq/data/alignment/yeast_genome/ ~/data/alignment/yeast_genome

cd ~/data/alignment/yeast_genome

# --simple means use subspecies strain name as name
perl         ~/Scripts/withncbi/taxon/strain_info.pl \
    --file   yeast_ncbi.csv                          \
    --simple \
    --id     559292                                  \
    --id     285006                                  \
    --id     307796                                  \
    --id     658763                                  \
    --id     574961                                  \
    --id     643680                                  \
    --id     721032                                  \
    --id     929585                                  \
    --id     226125                                  \
    --name   226125=Spar

# for paralog, need full chromosomes
perl       ~/Scripts/withncbi/taxon/strain_info.pl \
    --file yeast_for_paralog.csv                   \
    --id   658763                                  \
    --id   559292                                  \
    --id   285006                                  \
    --id   252598                                  \
    --id   1087981                                 \
    --id   889517                                  \
    --id   580240                                  \
    --id   28985                                   \
    --id   931890                                  \
    --id   284811                                  \
    --id   1071383                                 \
    --id   1071382                                 \
    --id   559295                                  \
    --id   226302                                  \
    --id   284593                                  \
    --id   1071378                                 \
    --id   1064592                                 \
    --id   1071380                                 \
    --id   1071381                                 \
    --id   4950                                    \
    --id   4956                                    \
    --id   1294331                                 \
    --id   1434269                                 \
    --name 658763=Sigma1278b                       \
    --name 559292=S288c                            \
    --name 285006=RM11_1a                          \
    --name 252598=Sbou                             \
    --name 1087981=YJSH1                           \
    --name 889517=CEN_PK113_7D                     \
    --name 580240=W303                             \
    --name 28985=Klac                              \
    --name 931890=Ecym                             \
    --name 284811=Agos                             \
    --name 1071383=Knag                            \
    --name 1071382=Kafr                            \
    --name 559295=Lthe                             \
    --name 226302=Lklu                             \
    --name 284593=Cgla                             \
    --name 1071378=Ndai                            \
    --name 1064592=Ncas                            \
    --name 1071380=Tbla                            \
    --name 1071381=Tpha                            \
    --name 4950=Tdel                               \
    --name 4956=Zrou                               \
    --name 1294331=YJM993                          \
    --name 1434269=UFMG_A_905

#----------------------------#
# yeast
#----------------------------#
cd ~/data/alignment/

# Saccharomyces cerevisiae
# copy sequences (raw)
perl           ~/Scripts/withncbi/taxon/strain_bz.pl        \
    --file     ~/data/alignment/yeast_genome/yeast_ncbi.csv \
    -w         ~/data/alignment                             \
    --seq_dir  ~/data/alignment/yeast_genome                \
    --name     yeast_new                                    \
    --use_name \
    -t         S288c                                        \
    -q         RM11_1a                                      \
    -q         YJM789                                       \
    -q         Sigma1278b                                   \
    -q         JAY291                                       \
    -q         EC1118                                       \
    -q         Kyokai_no_7                                  \
    -q         T7                                           \
    -q         Spar                                         \
    -o         Spar                                         \
    --keep     S288c                                        \
    --keep     YJM789                                       \
    --keep     Sigma1278b                                   \
    --keep     JAY291                                       \
    --keep     EC1118                                       \
    --keep     Kyokai_no_7                                  \
    --keep     T7                                           \
    --keep     Spar
    
sh yeast_new/1_real_chr.sh
sh yeast_new/2_file_rm.sh

# OPTIONAL
# copy RepeatMaskered folder back to ~/data/alignment/yeast_genome/

# don't copy sequences (RepeatMasker done)
perl           ~/Scripts/withncbi/taxon/strain_bz.pl        \
    --file     ~/data/alignment/yeast_genome/yeast_ncbi.csv \
    -w         ~/data/alignment                             \
    --name     yeast_new                                    \
    --use_name \
    -t         S288c                                        \
    -q         RM11_1a                                      \
    -q         Sigma1278b                                   \
    -q         YJM789                                       \
    -q         Spar                                         \
    -o         Spar


#----------------------------#
# yeast self
#----------------------------#
# The genome should be assembled at chromosome level.

cd ~/data/alignment/self

# Saccharomyces cerevisiae
# copy sequences (raw)
perl           ~/Scripts/withncbi/taxon/strain_bz_self.pl          \
    --file     ~/data/alignment/yeast_genome/yeast_for_paralog.csv \
    -w         ~/data/alignment/self                               \
    --seq_dir  ~/data/alignment/yeast_genome                       \
    --name     yeast_new                                           \
    --use_name \
    -t         S288c                                               \
    -q         Ecym                                                \
    -q         Agos                                                \
    -q         Knag                                                \
    -q         Kafr                                                \
    -q         Klac                                                \
    -q         Lthe                                                \
    -q         Lklu                                                \
    -q         Cgla                                                \
    -q         Ndai                                                \
    -q         Ncas                                                \
    -q         Sigma1278b                                          \
    -q         RM11_1a                                             \
    -q         Sbou                                                \
    -q         YJSH1                                               \
    -q         CEN_PK113_7D                                        \
    -q         W303                                                \
    -q         YJM993                                              \
    -q         UFMG_A_905                                          \
    -q         Tbla                                                \
    -q         Tpha                                                \
    -q         Tdel                                                \
    -q         Zrou

    
sh yeast_new/1_real_chr.sh
sh yeast_new/2_file_rm.sh

# OPTIONAL
# copy RepeatMaskered folder back to ~/data/alignment/yeast_genome/

# don't copy sequences (RepeatMasker done)
perl           ~/Scripts/withncbi/taxon/strain_bz_self.pl          \
    --file     ~/data/alignment/yeast_genome/yeast_for_paralog.csv \
    -w         ~/data/alignment/self                               \
    --name     yeast_new                                           \
    --use_name \
    -t         S288c                                               \
    -q         Ecym                                                \
    -q         Agos                                                \
    -q         Knag                                                \
    -q         Kafr                                                \
    -q         Klac                                                \
    -q         Lthe                                                \
    -q         Lklu                                                \
    -q         Cgla                                                \
    -q         Ndai                                                \
    -q         Ncas                                                \
    -q         Sigma1278b                                          \
    -q         RM11_1a                                             \
    -q         Sbou                                                \
    -q         YJSH1                                               \
    -q         CEN_PK113_7D                                        \
    -q         W303                                                \
    -q         YJM993                                              \
    -q         UFMG_A_905                                          \
    -q         Tbla                                                \
    -q         Tpha                                                \
    -q         Tdel                                                \
    -q         Zrou

perl           ~/Scripts/withncbi/taxon/strain_bz_self.pl          \
    --file     ~/data/alignment/yeast_genome/yeast_for_paralog.csv \
    -w         ~/data/alignment/self                               \
    --name     yeast_new                                           \
    --use_name \
    -t         S288c
