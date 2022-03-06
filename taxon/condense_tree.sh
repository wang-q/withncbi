#!/usr/bin/env bash

USAGE="
Usage: $0 NWK_FILE MAP_FILE NODE_COL GROUP_COL

Default values:
    NODE_COL    1
    GROUP_COL   2

$ bash condense_tree.sh tree.newick taxon.tsv
"

if [ "$#" -lt 2 ]; then
    echo >&2 "$USAGE"
    exit 1
fi

# Check whether newick-utils is installed
hash nw_clade 2>/dev/null || {
    echo >&2 "newick-utils is required but it's not installed.";
    echo >&2 "Install with homebrew: brew install brewsci/bio/newick-utils";
    exit 1;
}

# Set default parameters
NWK_FILE=$1
MAP_FILE=$2
NODE_COL=${3:-1}
GROUP_COL=${4:-2}

# Run
>&2 cat <<EOF
==> Parameters
NWK_FILE=${NWK_FILE}
MAP_FILE=${MAP_FILE}
NODE_COL=${NODE_COL}
GROUP_COL=${GROUP_COL}

Writing condense.newick and condense.map

EOF

mytmpdir=`mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir'`

nw_labels ${NWK_FILE} > ${mytmpdir}/labels.lst

cat $MAP_FILE |
    grep -v "^#" |
    tr " " "_" |
    tsv-join -k 1 -f ${mytmpdir}/labels.lst \
    > ${mytmpdir}/valid.tsv

cat ${mytmpdir}/valid.tsv |
    tsv-select -f $GROUP_COL |
    tsv-uniq |
    grep -v "NA" |
    grep "." \
    > ${mytmpdir}/group.lst

# Check monophyly for group
for GROUP in $(cat ${mytmpdir}/group.lst); do
    cat ${mytmpdir}/valid.tsv |
        tsv-filter --str-eq "$GROUP_COL:$GROUP" |
        tsv-select -f $NODE_COL \
        > ${mytmpdir}/${GROUP}.lst

    NODE=$(
        nw_clade -m ${NWK_FILE} $(cat ${mytmpdir}/${GROUP}.lst) |
            nw_stats -f l - |
            cut -f 3
    )

    if [[ "$NODE" ]]; then
        echo "${GROUP}" >> ${mytmpdir}/group.monophyly.lst
        cat ${mytmpdir}/${GROUP}.lst |
            xargs -I[] echo "[] ${GROUP}___${NODE}" \
            >> ${mytmpdir}/group.monophyly.map
    fi

done

# Merge labels in group to a higher-rank
nw_rename ${NWK_FILE} ${mytmpdir}/group.monophyly.map |
    nw_condense - \
    > ${mytmpdir}/group.monophyly.newick

# Outputs
mv ${mytmpdir}/group.monophyly.newick condense.newick
mv ${mytmpdir}/group.monophyly.map condense.map

rm -fr ${mytmpdir}
