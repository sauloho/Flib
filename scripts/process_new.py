#!/usr/bin/env python
#
# Take a Rosetta++ fragment library and generate a valid SAINT2
# fragment library.
# 
# Format:
#
# F <frag_id> P <pos> L <length> S <score> = <info>
# <res_id> <phi> <psi> <omega> <N angle> <CA angle> <C angle>
# <res_id> <phi> <psi> <omega> <N angle> <CA angle> <C angle>
# <res_id> <phi> <psi> <omega> <N angle> <CA angle> <C angle>
#
# Example:
# 
#
# F 0 P 0 L 9 S 1.511 = 1fcd A 116 R 0
# 0  -77.7  143.6  179.6  126.7  115.1  114.9 
# 1 -145.6  173.3  173.7  123.0  107.1  111.6 
# 2  -57.2  -38.7 -179.5  122.2  107.5  115.7 
# 3  -86.8   -8.1  175.3  122.6  113.2  118.0 
# 4  -82.5  -26.7  176.9  120.5  107.8  115.5 
# 5  -44.1  -39.7 -178.4  126.1  113.1  117.2 
# 6  -68.7  -14.9 -176.2  123.9  118.0  116.8 
# 7 -119.8  -34.2  177.6  123.8  112.5  114.4 
# 8 -124.7   95.4 -171.4  122.6   99.9  121.6 
#
import math
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Vector import calc_dihedral, calc_angle
from sys import stdin, stderr, argv, exit
from angles import print_pp

def make_filename(pdb):
    ### Location of PDB database: ###
    pdb=pdb.lower()
    return argv[1]+"/"+ pdb+".pdb"

def get_pp(pdb, chain,start, length,seq ):
	"""retrieve the residiues for a given pdb file and chain as polypeptides"""
	f=make_filename(pdb)
	p=PDBParser(PERMISSIVE=1)
	pdb_struct = p.get_structure(pdb, f) 		# Load the pdb structure pdb contained on the file f.
	pdb_chain = pdb_struct[0][chain]	  		# Select the right Chain of the structure.	
	ppb = PPBuilder()					  		# Initialize a peptide builder.
	peptides = ppb.build_peptides(pdb_chain)	# Load the given chain as a peptide.	
	for i,pep in enumerate(peptides):
		if str(pep.get_sequence()).find(seq) != -1 :
			start = str(pep.get_sequence()).find(seq)
			break
	if start > 0 and (start+length+2) <= len(pep):
		pp=pep[(start-1):(start+length+2)]
		return pp
	else:
		raise

fails = 0
fails2 = 0
s = 0 # Fragment "serial number".

if __name__ == '__main__':
	for line in stdin:
		frag = line.split("\t")
# 0		  1		   2	   3	  4				  5		  6	      7		  8			9	  10 	  11	  12
# 3DXX    A        65      73     VGACLASAH       H       0       1       9         0     2.050   10      1.780
		pdb   = frag[0].lower()		
		chain = frag[1]
		start = int(frag[2])
		end = frag[3]
		seq = frag[4]
		length = int(frag[8])
		pos = frag[9]
		score = "%.2f" % ( math.exp((float(frag[12].rstrip()))/1000))
	  	if float(score)<0.1 or math.isnan(float(score)):
			score="1.0"
            
		try:
			pp = get_pp(pdb, chain,start,length,seq) 			# Polypeptide pp now contains atomic information for fragment.
		except:
			print >> stderr, "E: failed to process", pdb, chain, start # Cannot open file.
			fails += 1
		else:
			str1=""
			for res in pp[1:length+1]:			
				 str1+=three_to_one(res.get_resname())
			if str(seq[0:length]) == str(str1):
				#      F   0   P    0    L   9        S   1.511    =   1fcd    A 116 R 0
				print "F ",s," P ",pos," L ",str(length)," S ",score," = ",pdb," ",chain,"    ",str(start)#
				try:
					print_pp(pp,offset=1,lng=length)
					s+=1
					if length > 6:
						print "F ",s," P ",pos," L ",str(6)," S ",score," = ",pdb," ",chain,"    ",str(start)
						print_pp(pp,offset=1,lng=6)
						s+=1
				except:
					continue
			elif not seq is str1:
				fails2+=1
				#print >> stderr,  seq+" "+str1
				print >> stderr, "E: Sequence mismatch:", pdb, chain, start
print >> stderr, "%d frags extracted from PDB" % s
print >> stderr, "%d processing failures" % fails
print >> stderr, "%d sequence failures" % fails2

