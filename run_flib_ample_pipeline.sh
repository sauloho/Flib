#!/bin/bash

# ITERATE OVER A RANGE OF LENGTH INTERVALS OF SIZE 5 AND PRODUCE FRAGMENTS USING THIS LENGTH SETTING:
for I in {5..45..5}
do
    # FOR A DETAILED LIST OF THE COMMAND OPTIONS FOR FLIB, USE OPTION --help
    ./src/Flib --coevo_only -i $1 -C $1.con -M 200 -l $I -L $[I+5] 2> /dev/null > $1.flib_ample_"$I"
    cp  $1.flib_ample_"$I" $1.flib_ample_"$I"_nh
    for HOMOLOG in $(cat $1.homol)
    do
        sed -e "/$HOMOLOG/d" "$1".flib_ample_"$I"_nh > "$1".tmp
        mv $1.tmp $1.flib_ample_"$I"_nh
    done 
done

cat $1.flib_ample_*_nh | sort -k 13,13n > $1.lib_reduced
