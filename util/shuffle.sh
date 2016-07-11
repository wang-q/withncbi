#!/usr/bin/env bash

USAGE="Usage: $0 FASTA_FILE SEQ_NUMBER MINIMAL_LENGTH MAXIMAL_LENGTH"

# bash shuffle.sh ~/data/rna-seq/medfood/process/Arabidopsis_thaliana/result/Trinity.fasta 100 1000 2000

if [ "$#" -lt 1 ]; then
    echo >&2 "$USAGE"
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
SEQ_NUMBER=${2:-100}
MINIMAL_LENGTH=${3:-1000}
MAXIMAL_LENGTH=${4:-2000}

echo >&2 "==> Parameters <=="
echo >&2 "    FASTA_FILE=${FASTA_FILE}"
echo >&2 "    SEQ_NUMBER=${SEQ_NUMBER}"
echo >&2 "    MINIMAL_LENGTH=${MINIMAL_LENGTH}"
echo >&2 "    MAXIMAL_LENGTH=${MAXIMAL_LENGTH}"

echo >&2 "==> Getting list of sequences"
faops filter -a ${MINIMAL_LENGTH} -z ${MAXIMAL_LENGTH} ${FASTA_FILE} stdout \
    | perl -nl -e '/^>(\S+)/ or next; print $1' \
    | perl -MList::Util=shuffle -wne 'print shuffle <>;' \
    | head -n ${SEQ_NUMBER} \
    | sort \
    > seq_list.tmp

[ $? -ne 0 ] && echo >&2 "Getting list of sequences failed" && exit 255

echo >&2 "==> Getting sequences from list"
faops some ${FASTA_FILE} seq_list.tmp stdout
