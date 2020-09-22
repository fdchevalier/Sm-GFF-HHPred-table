#!/bin/bash
# Title: hhpred-ann.sh
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2015-05-23
# Modified in: 2020-09-21
# Licence : GPL v3



#======#
# Aims #
#======#

aim="Run HHPred on a list of sequences and collect best annotation in a table."



#==========#
# Versions #
#==========#

# v1.0 - 2020-09-21: reorganization
# v0.0 - 2015-05-23: creation

version=$(grep -i -m 1 "version" "$0" | cut -d ":" -f 2 | sed "s/^ *//g")



#===========#
# Functions #
#===========#

# Usage message
function usage {
    echo -e "
    \e[32m ${0##*/} \e[00m -i|--in file -d|--db path -o|--out directory -p|--plim integer -h|--help

Aim: $aim

Version: $version

Options:
    -i, --in        input list file
    -d, --db        path to database. Several databases can be specified at once.
    -o, --out       output directory (optional) [default: results]
    -p, --plim      probability threshold value (optional) [default: 50]
    -h, --help      this message
    "
}


# Info message
function info {
    if [[ -t 1 ]]
    then
        echo -e "\e[32mInfo:\e[00m $1"
    else
        echo -e "Info: $1"
    fi
}


# Warning message
function warning {
    if [[ -t 1 ]]
    then
        echo -e "\e[33mWarning:\e[00m $1"
    else
        echo -e "Warning: $1"
    fi
}


# Error message
## usage: error "message" exit_code
## exit code optional (no exit allowing downstream steps)
function error {
    if [[ -t 1 ]]
    then
        echo -e "\e[31mError:\e[00m $1"
    else
        echo -e "Error: $1"
    fi

    if [[ -n $2 ]]
    then
        exit $2
    fi
}


# Dependency test
function test_dep {
    which $1 &> /dev/null
    if [[ $? != 0 ]]
    then
        error "Command $1 not found. Exiting..." 1
    fi
}


# Progress bar
## Usage: ProgressBar $mystep $myend
function ProgressBar {
    if [[ -t 1 ]]
    then
        # Process data
        let _progress=(${1}*100/${2}*100)/100
        let _done=(${_progress}*4)/10
        let _left=40-$_done
        # Build progressbar string lengths
        _fill=$(printf "%${_done}s")
        _empty=$(printf "%${_left}s")

        # Build progressbar strings and print the ProgressBar line
        # Output example:
        # Progress : [########################################] 100%
        #printf "\rProgress : [${_fill// /=}${_empty// / }] ${_progress}%%"
        printf "\r\e[32mProgress:\e[00m [${_fill// /=}${_empty// / }] ${_progress}%%"
        
        [[ ${_progress} == 100 ]] && echo ""
    fi
}


# Clean up function for trap command
## Usage: clean_up file1 file2 ...
function clean_up {
    rm -rf $@
    exit 1
}



#==============#
# Dependencies #
#==============#

test_dep hhsearch



#===========#
# Variables #
#===========#

# Options
while [[ $# -gt 0 ]]
do
    case $1 in
        -i|--in     ) list="$2" ; shift 2 ;;
        -d|--db     ) mydb="-d $(printf '%q' "$2")" ; shift 2 
                        while [[ ! -z "$1" && $(echo "$1"\ | grep -qv "^-" ; echo $?) == 0 ]]
                        do
                            mydb="$mydb -d $(printf '%q' "$1")"
                            shift
                        done ;; 
        -o|--out    ) output="$2" ; shift 2 ;;
        -p|--plim   ) plim="$2" ; shift 2 ;;
        -s|--switch ) switch=1 ; shift ;;
        -h|--help   ) usage ; exit 0 ;;
        *           ) error "Invalid option: $1\n$(usage)" 1 ;;
    esac
done


# Check the existence of obligatory options
[[ -z "$list" ]] && error "The option input is required. Exiting...\n$(usage)" 1

# Set defaults for optional options
[[ -z "$outdir" ]] && outdir="results"
[[ -z "$plim" ]]   && plim=50
# elim="$3"

# Output filename
annfile="${list##*/}.db"



#============#
# Processing #
#============#

# Check if output file exists
[[ -e "$outdir/$annfile" ]] && echo "$outdir/$annfile exists already. Exiting..." && exit 1

# Create output directory if needed
[[ ! -d "$outdir" ]] && mkdir -p "$outdir"

# Analyze each sequence from the list
count=0
wc_ls=$(wc -l < "$list")
while read i
do
    # Counter
    ((count++))
    ProgressBar $count $wc_ls

    # Transcript and gene names
    name="${i##*/}"
    trans=$(echo "${name%.*}" | sed s"/_/./2")

    # HHPred analysis
    hhsearch -cpu $(nproc) -i "$i" $mydb -o "$outdir/${name%.*}.hhr" -v 0 
    
    # Default annotation
    ann="NA"

    # Probability column values
    # eval=$(sed -n "/No Hit/,/No 1/{//d;p}" "$outdir/${name%.*}.hhr" | sed "/^$/d" | sed -r "s/.{41}(.{8}).*/\1/")
    pval=$(sed -n "/No Hit/,/No 1/{//d;p}" "$outdir/${name%.*}.hhr" | sed "/^$/d" | sed -r "s/.{35}(.{5}).*/\1/")

    # List of annotations from alignment entries
    list_ann=$(grep ">" "$outdir/${name%.*}.hhr" | sed "s/>[a-zA-Z0-9_.]* //")
    for ((i=1; i <= $(wc -l <<<"$list_ann"); i++))
    do
        # Test probability value
        # [[ $(awk -v a="$(sed -n "${i}p" <<<"$eval")" -v elim="$elim" 'BEGIN{print (a < elim)?1:0}') -eq 0 ]] && continue
        [[ $(awk -v a="$(sed -n "${i}p" <<<"$pval")" -v plim="$plim" 'BEGIN{print (a >= plim)?1:0}') -eq 0 ]] && continue

        # Get annotation
        ann_tmp=$(sed -n "${i}p" <<<"$list_ann" | cut -d ";" -f 1)

        # Check annotation
        ## Is there any letter?
        if [[ ! $(grep "[[:alpha:]]" <<<"$ann_tmp") ]]
        then
            NF=$(sed -n "${i}p" <<<"$list_ann" | awk -F ";" 'END {print NF}')
            for ((j=2; j<=$NF; j++))
            do
                ann_tmp=$(sed -n "${i}p" <<<"$list_ann" | cut -d ";" -f $j)
                [[ $(grep "[[:alpha:]]" <<<"$ann_tmp") ]] && break
            done
        fi

        # Is annotation a sequence?
        if [[ ! $(egrep "[A-Z]{60}" <<<"$ann_tmp") ]]
        then
            # Remove unnecessary space 
            ann=$(sed "s/^ *//" <<<"$ann_tmp")

            # Clean SCOPe annotation
            ann=$(sed -r "s/^[a-z]\.[0-9]*\.[0-9]*\.[0-9]* \([^\)]*\) //g" <<<"$ann")

            break
        fi
    done

    echo -e "$trans\t$ann" >> "$outdir/$annfile"
done < "$list"

exit 0
