#!/bin/bash

# SET UP THE PATHS FOR REQUIRED DEPENDENCIES:
export spineXblast=/data/cockatrice/wilman/progs/ncbi-blast-2.2.27+/
export spineXcodir=/users/oliveira/SPINEX/spineXpublic/
export PSIPRED=/users/oliveira/PSIPRED
export BLAST=/users/oliveira/blast/blast-2.2.17/bin/
export BLASTDB=/users/oliveira/PSIPRED/unirefdb90
export BLAST_PDB=/users/oliveira/blast/blast-2.2.17/db/pdb_seqres.txt
export HHSUITE=/users/oliveira/hhsuite-2.0.16-linux-x86_64/
export HHBLITSDB=$HHSUITE/hhblits_database/pdb70_05Jun14
export FLIB=/users/oliveira/Flib/
export PDB=/users/oliveira/PDB/

# READ ARGUMENTS:
OUTPUT=$1
#read OUTPUT

# ADJUST THE FILES THAT WILL BE GENERATED/COMPUTED:
generate_ss=false		# Change value to "true" if running local version of PSIPRED.
generate_pssm=false		# Change value to "true" if running local version of SPINE-X.
generate_spinex=false	# Change value to "true" if running local version of SPINE-X.
generate_hhr=true		# Change value to "true" if running local version of HHBlits.	
generate_flib=false		# Change value to "true" if running local version of Flib. 
parse_flib=false 		# Change value to "true" if parsing Flib libraries to SAINT2 format. 
remove_homologs=false	# Change value to "true" if removing homologs from frag. libraries.

##### GENERAL SET UP #####

# Generate .ss file
if [ "$generate_ss" = true ] ; then
	echo "--------------------------------------"
	echo "Generating Secondary Structure Prediction file:"
	$PSIPRED/runpsipredplus ./$OUTPUT.fasta.txt
	echo "Done"
	echo "--------------------------------------"
fi

# Generate .mat file
if [ "$generate_pssm" = true ] ; then
	echo "Generating PSSM file for SPINE-X torsion angle prediction:"
	$spineXblast/bin/psiblast -db $BLASTDB -out_pssm $OUTPUT.txt -evalue 0.02 -query ./$OUTPUT.fasta.txt -out_ascii_pssm ./$OUTPUT.mat -out $OUTPUT.out -num_iterations 5
	echo "Done"
	echo "-------------------------------------"
fi

# Generate SPINE-X file
if [ "$generate_spinex" = true ] ; then
	echo "Generating SPINE-X torsion angle prediction file:"
	echo $OUTPUT > $OUTPUT.temp
	$spineXcodir/spX.pl ./$OUTPUT.temp ./
	mv ./spXout/$OUTPUT.spXout ./
	rm $OUTPUT.temp
	echo "Done"
	echo "--------------------------------------"
fi

# Run HHSearch
if [ "$generate_hhr" = true ] ; then
        echo "Generating HHSearch File:"
        $HHSUITE/bin/hhblits -d $HHBLITSDB -i ./$OUTPUT.fasta.txt -o $OUTPUT.hhr
        echo "Done"
fi

# Generate list of Homologs 
if [ "$remove_homologs" = true ] ; then
       # Generate list of homologs
       $BLAST/blastall -p blastp -i ./$OUTPUT.fasta.txt -d $BLAST_PDB -e 0.05 -m 8  > $OUTPUT.blast
       cat $OUTPUT.blast | awk '{print substr($2,1,4)}' | sort | uniq > $OUTPUT.homol 
fi

##### FLIB #####
if [ "$generate_flib" = true ] ; then
	echo "Generating FLIB File:"
   	$FLIB/Flib $OUTPUT $PDB > $OUTPUT.lib3000 2> $OUTPUT.log           # Generates LIB3000
	sort -k 10,10n -k 13,13n $OUTPUT.lib3000 > $OUTPUT.tmp;                   # Sorts LIB3000
	mv $OUTPUT.tmp $OUTPUT.lib3000;                                           

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
	awk '$6 == "H" {print }' $OUTPUT.lib_final > $OUTPUT.combined
    awk '$6 == "B" {print }' $OUTPUT.lib20_ori >> $OUTPUT.combined
    awk '$6 == "O" {print }' $OUTPUT.lib_final >> $OUTPUT.combined
    awk '$6 == "O" {print }' $OUTPUT.lib20	   >> $OUTPUT.combined
    awk '$6 == "L" {print }' $OUTPUT.lib_final >> $OUTPUT.combined
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



