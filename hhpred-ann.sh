#!/bin/bash
# Created in: 2015-05-23


#===========#
# Variables #
#===========#

list="$1"
clist="list"
outdir="results"

annfile="${list##*/}.db"

# Environment
export HHLIB="$HOME/local/bin/hhsuite/"
export HHDB="$HOME/bio/HH-suite/hhsearch_dbs/"

if [[ ! -f "$clist" ]]
then
    echo "Complete list unavailable"
    exit 1
fi

if [[ ! -d "$outdir" ]]
then
    mkdir "$outdir"
fi


#============#
# Processing #
#============#

for i in $(cat $list)
do
    name="${i##*/}"
    gene=$(echo "$name" | cut -d "_" -f -2)

    # Output gene name
    if [[ $(grep -c "$gene" "$clist") -eq 1 ]]
    then
        genenb="$gene"
    else
        genenb=$(echo "${name%.*}" | sed s"/_/./2")
    fi

    # Analysis
    hhsearch -cpu $(nproc) -i "$i" -d "$HHDB/pdb70_06Sep14.hhm" -o "$outdir/${name%.*}.hhr" 1> /dev/null
    
    # Annotation
    ann=$(grep -m 1 ">" "$outdir/${name%.*}.hhr" | sed "s/>[a-zA-Z0-9_]* //;s/;.*//")
    echo -e "$genenb\t$ann" >> "$outdir/$annfile"
done

exit 0
