#!/usr/bin/env python

import os,sys
import re

if __name__ == '__main__':
	counts = open(sys.argv[1]+".hhr",'r')
#	homols = open(sys.argv[1]+".homol",'r')
	homo_list = []
#	for homolog in homols:
#		homo_list.append(homolog.rstrip())
#	homols.close()
	#print homo_list
	i=1000
	for line in counts:
		if line[0] == '>':
			i=0
		if i==1:
			m=re.search("Score=([0-9]+.[0-9]+)",line)
			score = m.group(1)					
		if i==8:
			aux = line.split()
			ss_seq1 = aux[2].upper()			
		if i==3:
			line = line[2:]
			line = line[0:4] + "          " + line[14:]
			seq1=line
		if i==7:
			line = line[2:]
			pdb = line[0] + line[1] + line[2] + line[3]
			line = line[0:4] + " " + line[5:]
			seq2 = line
		if i==9:
			aux = line.split()
			ss_seq2 = aux[2].upper()			
			if homo_list.count(pdb) == 0:
				seq1_list = seq1.split()
				seq1_list[1] = int(seq1_list[1])
				seq2_list = seq2.split()
				seq2_list[2] = int(seq2_list[2])
				for cont1 in range(0,len(seq1_list[2])-9):
					begin1 = seq1_list[1] - 1 + cont1
					begin2 = seq2_list[2] - 1 + cont1
					output=seq2_list[0]+'\t'+seq2_list[1]+'\t'+str(begin2)+'\t'
					if seq1_list[2][cont1:cont1+9].count("-")==0 and seq2_list[3][cont1:cont1+9].count("-")==0:
						output=output+str(seq2_list[2]-1+cont1+8)+'\t'+ seq2_list[3][cont1:cont1+9] +"\t"
						output=output+'X'+"\t0\t0\t9\t"+str(begin1)+"\t"+score+"\t0"
						output=output+"\t"+seq1_list[2][cont1:(cont1+9)]+"\t"+ss_seq1[cont1:(cont1+9)]+"\t"+ss_seq2[cont1:(cont1+9)]
						print output
#Header,&Chain,&start1,&k,DB_Seq,&c,&ramach_score,&seq_score,&length,&start2,&resolution,&ss_score
		i=i+1
