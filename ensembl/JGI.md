# JGI data

## Scripts

Same for each species.

### `jgi_fasta.sh`

Genome, RepeatMasker and gff3.

```bash

cat <<'EOF' > ~/data/alignment/JGI/jgi_fasta.sh
#!/bin/bash

USAGE="Usage: $0 GENOME_NAME"

if [ "$#" -lt 1 ]; then
    echo "$USAGE"
    exit 1
fi

echo "==> parameters <=="
echo "    " $@

GENOME_NAME=$1

if [ ! -d ~/data/alignment/JGI/${GENOME_NAME} ]
then
    echo "==> ${GENOME_NAME}"

    mkdir -p ~/data/alignment/JGI/${GENOME_NAME}
    cd ~/data/alignment/JGI/${GENOME_NAME}
    
    find ~/data/PhytozomeV11/${GENOME_NAME}/assembly -type f -name "*.fa.gz" -not -name "*masked" \
        | xargs gzip -d -c > toplevel.fa
    
    faops count toplevel.fa \
        | perl -aln -e '
            next if $F[0] eq 'total';
            print $F[0] if $F[1] > 50000;
            print $F[0] if $F[1] > 5000  and $F[6]/$F[1] < 0.05;
            ' \
        | uniq > listFile
    faops some toplevel.fa listFile toplevel.filtered.fa
    faops split-name toplevel.filtered.fa .
    rm toplevel.fa toplevel.filtered.fa listFile

    cat scaffold*.fa > Un.fa
    rm scaffold*.fa
    mv Un.fa Un.fa.skip
else
    echo "==> ~/data/alignment/JGI/${GENOME_NAME} exists"
fi

EOF

```

## Athaliana 

```bash
mkdir -p ~/data/alignment/JGI
cd ~/data/alignment/JGI

bash ~/data/alignment/JGI/jgi_fasta.sh Athaliana
```

