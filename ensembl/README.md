# Ensembl related scripts.

Configurations stored in `ensembl_82.yml`.

```bash
cd ~/Scripts/withncbi/ensembl/
perl ensembl_batch.pl -i ensembl_82.yml

bash ensembl.build.sh
bash ensembl.fasta.sh
bash ensembl.anno.sh

cp ensembl.initrc.pm ~/Scripts/alignDB/
```
