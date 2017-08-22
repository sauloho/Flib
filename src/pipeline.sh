#!/bin/bash

# SET UP THE PATHS FOR REQUIRED DEPENDENCIES:
export PSIPRED=/home/U055415/PSIPRED
export SPIDER=/home/U055415/SPIDER2_local/misc/
export HHSUITE=/home/U055415/hh-suite/bin/
export HHBLITSDB=$HHSUITE/../db/pdb70
export BLAST=/usr/bin/
export BLAST_PDB=/home/U055415/Databases/pdbaa
export FLIB=/home/U055415/Flib/
export PDB=/run/media/U055415/DATADRIVE0/PDB/

# READ ARGUMENTS:
OUTPUT=$1

# ADJUST THE FILES THAT WILL BE GENERATED/COMPUTED:
generate_ss=true		# Change value to "true" if running local version of PSIPRED.
generate_spider2=true		# Change value to "true" if running local version of SPINE-X.
generate_hhr=true		# Change value to "true" if running local version of HHBlits.	
generate_flib=true		# Change value to "true" if running local version of Flib. 
parse_flib=true 		# Change value to "true" if parsing Flib libraries to SAINT2 format. 
remove_homologs=true		# Change value to "true" if removing homologs from frag. libraries.

##### GENERAL SET UP #####

# Generate SS Prediction using PSIPRED
if [ "$generate_ss" = true ] ; then
	echo "--------------------------------------"
	echo "Generating Secondary Structure Prediction using PSIPRED:"
	$PSIPRED/runpsipredplus ./$OUTPUT.fasta.txt
	echo "Done"
	echo "--------------------------------------"
fi

# Generate SPIDER2 Torsion Angle Prediction
if [ "$generate_spider2" = true ] ; then
	echo "Generating SPIDER2 Torsion Angle Prediction:"
	cat $OUTPUT.fasta.txt > $OUTPUT.seq
	$SPIDER/run_local.sh $OUTPUT.seq 2> $OUTPUT.spd_err
	rm $OUTPUT.seq
	echo "Done"
	echo "--------------------------------------"
fi

# Run HHSearch
if [ "$generate_hhr" = true ] ; then
        echo "Generating HHSearch File:"
        $HHSUITE/hhblits -d $HHBLITSDB -i ./$OUTPUT.fasta.txt -o $OUTPUT.hhr
        echo "Done"
fi

# Generate list of Homologs 
if [ "$remove_homologs" = true ] ; then
       # Generate list of homologs
       blastp -query ./$OUTPUT.fasta.txt -db $BLAST_PDB -evalue 0.05 -outfmt 6  > $OUTPUT.blast
       awk '{print $2}' $OUTPUT.blast | sed -e "s/.*pdb|//g" | cut -c 1-4 | sort | uniq | tr '[:upper:]' '[:lower:]' > $OUTPUT.homol
fi

##### FLIB #####
if [ "$generate_flib" = true ] ; then
	echo "Generating FLIB File:"
   	$FLIB/Flib $OUTPUT > $OUTPUT.lib3000 2> $OUTPUT.log           # Generates LIB3000
	sort -k 10,10n -k 13,13n $OUTPUT.lib3000 > $OUTPUT.tmp;                   # Sorts LIB3000
	mv $OUTPUT.tmp $OUTPUT.lib3000;                                           
	cp $OUTPUT.lib3000 $OUTPUT.lib3000_ori

	### HOMOLOG REMOVAL ####
	if [ "$remove_homologs" = true ] ; then
		cp $OUTPUT.lib3000 $OUTPUT.lib3000_nh;
		for HOMOLOG in $(cat $OUTPUT.homol)
		do
			sed -e "/$HOMOLOG/d" "$OUTPUT".lib3000_nh > "$OUTPUT".tmp
		        mv $OUTPUT.tmp $OUTPUT.lib3000_nh
		done
		mv $OUTPUT.lib3000_nh $OUTPUT.lib3000
	fi

	$FLIB/filterlib2 $OUTPUT $OUTPUT.lib3000                             # This will generate LIB20 and LIB500
	mv $OUTPUT.lib20 $OUTPUT.lib20_ori
	mv $OUTPUT.lib500 $OUTPUT.lib500_ori
	
	### Fragment Library Enrichment ###
	$FLIB/Flib_Enrich $OUTPUT.lib500_ori $PDB 0.5 0 > $OUTPUT.clib 2> $OUTPUT.error # This will enrich LIB20 with fragments from LIB500: CLIB
	cat $OUTPUT.lib20_ori >> $OUTPUT.lib_tmp2                                     # Merges LIB20 and CLIB
	cat $OUTPUT.clib >> $OUTPUT.lib_tmp2                                      # ...
	sort -k 10,10n $OUTPUT.lib_tmp2 > $OUTPUT.lib_final                       # Sorts LIB20+CLIB : LIB_FINAL
	rm $OUTPUT.lib_tmp2     

        ### Parsing fragments from Threading hits: ###
	python $FLIB/parse_hhr.py $OUTPUT > $OUTPUT.lib_hhr 2> $OUTPUT.log    # Creates lib from threading hits: LIB_HHR
	sort -k 10,10n -k 11,11nr $OUTPUT.lib_hhr > $OUTPUT.ordered;                    # Sorts LIB_HHR
	$FLIB/parse_hhr $OUTPUT.ordered > $OUTPUT.lib9               # Parses LIB_HHR
	cut -f 1-12 $OUTPUT.lib9 > $OUTPUT.lib9_aux
	mv $OUTPUT.lib9_aux $OUTPUT.lib9
	awk 'BEGIN {OFS = "\t"} ; {print $0,"0.0"}' $OUTPUT.lib9 > $OUTPUT.rmsd_9_lib

        ### HOMOLOG REMOVAL ####
	if [ "$remove_homologs" = true ] ; then
		cp $OUTPUT.rmsd_9_lib $OUTPUT.rmsd_9_lib_nh;
	        for HOMOLOG in $(cat $OUTPUT.homol)
	        do
	       	        sed -e "/$HOMOLOG/d" "$OUTPUT".rmsd_9_lib_nh > "$OUTPUT".tmp
	              	mv $OUTPUT.tmp $OUTPUT.rmsd_9_lib_nh
	        done
	        mv $OUTPUT.rmsd_9_lib_nh $OUTPUT.rmsd_9_lib
	fi
	$FLIB/filterlib2 $OUTPUT $OUTPUT.rmsd_9_lib                 # Filter LIB_HHR so it does not contain more than 20 frags. per position.
	cat $OUTPUT.lib_final > $OUTPUT.combined
	awk '$6 == "O" {print }' $OUTPUT.lib20	   >> $OUTPUT.combined
	awk '$6 == "L" {print }' $OUTPUT.lib20     >> $OUTPUT.combined

	sort -k 10,10n -k 13,13nr $OUTPUT.combined > $OUTPUT.lib                        # The resulting library is LIB
	echo "Done"
	echo "--------------------------------------"
fi


##### LIBRARY PARSING #####

# Parse FLIB library into SAINT2 compliant format.
if [ "$parse_flib" = true ] ; then
	echo "Parsing Flib File:"
	python $FLIB/process_new.py $PDB < $OUTPUT.lib > $OUTPUT.flib 2> $OUTPUT.log
	echo "Done"
	echo "--------------------------------------"
fi



