#!/bin/bash

USAGE="Usage: $0 GENOME_NAME"

if [[ "$#" -lt 1 ]]; then
    echo "$USAGE"
    exit 1
fi

echo "==> parameters <=="
echo "    " $@

GENOME_NAME=$1
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

echo "==> mite_stat"
faops size mite.fa \
    | perl -nla -F"\t" -MStatistics::Descriptive -e '
        BEGIN {$stat = Statistics::Descriptive::Full->new;}
        next unless defined $F[1];
        $stat->add_data($F[1]);
        END {
            printf qq{Total:\t%d\nSum:\t%d\n}, $stat->count, $stat->sum;
            printf qq{Median:\t%d\nMean:\t%.1f\n}, $stat->median, $stat->mean;
            printf qq{Min:\t%d\nMax:\t%d\n}, $stat->min, $stat->max;
            printf qq{\n};
            printf qq{Length\t#Seqs\n};
            %distrib = $stat->frequency_distribution(10);
            for (sort {$a <=> $b} keys %distrib) {
                printf qq{%.1f\t%d\n}, $_, $distrib{$_};
            }
            printf qq{\n};
        }' \
    > mite_stat.txt

echo "==> genome blast"
faops filter -a 40  mite.fa mite.filter.fa
perl ~/Scripts/egaz/fasta_blastn.pl -f mite.filter.fa -g genome.fa -o mite.bg.blast --parallel 8
perl ~/Scripts/egaz/blastn_genome.pl -f mite.bg.blast -g genome.fa -o mite.bg.fasta -c 0.95 --parallel 8
cat mite.fa mite.bg.fasta \
    | faops filter -u stdin stdout \
    | faops filter -a 40 stdin stdout \
    > mite.all.fasta

echo "==> sparsemem_exact"
perl ~/Scripts/egaz/sparsemem_exact.pl -f mite.all.fasta -l 40 -g genome.fa -o mite.replace.tsv
cat mite.replace.tsv \
    | perl -nla -F"\t" -e ' print for @F' \
    | grep ':' \
    | sort | uniq \
    > mite.position.txt

echo "==> mite covers"
wc -l mite.position.txt >> mite_stat.txt
runlist cover mite.position.txt -o mite.yml

runlist span --op excise -n 10 mite.yml   -o mite.1.yml # remove small spans
runlist span --op fill   -n 10 mite.1.yml -o mite.2.yml # fill small holes
mv mite.2.yml mite.yml
rm mite.1.yml

runlist stat mite.yml -s chr.sizes

echo "==> clean"
rm mite.all.fasta mite.bg.blast mite.bg.fasta mite.filter.fa mite.replace.tsv

mv mite.yml ../yml
mv mite.yml.csv ../stat
mv mite_stat.txt ../stat
echo mite >> ../yml/repeat.family.txt
