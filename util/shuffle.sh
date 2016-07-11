#!/usr/bin/env bash

USAGE="Usage: $0 FASTA_FILE MINIMAL_LENGTH SEQ_NUMBER"

# bash shuffle.sh ~/data/rna-seq/medfood/process/Arabidopsis_thaliana/result/Trinity.fasta 1000 100

if [ "$#" -lt 1 ]; then
    echo "$USAGE"
    exit 1
fi

# chech faops is installed
hash faops 2>/dev/null || {
    echo >&2 "faops is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install wang-q/tap/faops";
    exit 1;
}

# set default parameters
FASTA_FILE=$1
MINIMAL_LENGTH=${2:-1000}
SEQ_NUMBER=${3:-100}

echo "==> Parameters <=="
echo "    FASTA_FILE=${FASTA_FILE}"
echo "    MINIMAL_LENGTH=${MINIMAL_LENGTH}"
echo "    SEQ_NUMBER=${SEQ_NUMBER}"

echo "==> Getting list of sequences"
faops filter -a ${MINIMAL_LENGTH} ${FASTA_FILE} stdout \
    | perl -nl -e '/^>(\S+)/ or next; print $1' \
    | perl -MList::Util=shuffle -wne 'print shuffle <>;' \
    | head -n ${SEQ_NUMBER} \
    | sort \
    > seq_list.tmp

[ $? -ne 0 ] && echo >&2 "Getting list of sequences failed" && exit 255

echo "==> Getting sequences from list"
faops some ${FASTA_FILE} seq_list.tmp stdout
