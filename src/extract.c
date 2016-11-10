/*****************************************************************/
/* 																 */		
/*     C. Deane - Superposition Function - 27/10/04				 */
/* 																 */
/*	Take in lists of fragments perform pariwise superpositon.	 */
/*	then cluster on the basis of these distances into families.  */
/* 																 */
/*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "slink.h"

#define MaxFrag 10000
#define TRUE 1
#define FALSE 0

#define min(a,b) (((a)<(b))?(a):(b))

BITS Atmrec;

/* Collect the coordinates from the pdb file starting from the correct residue. */
int extract_frags( char bum[FILEL],int length, int start_res,char ChainQ,char* Seq);

int main (int argc,char *argv[])
{
	/* this is a list of file names*/
	Data LLL;

	pid_t pid;
	int i,j,k,l,status,temp_nres,temp_natom;
	int h,b,o,top_max,cont;
	int length,index,numfrags=0;
	int start1;
	char Header[82],Chain, aux;
	char AUX[82],AUX2[82];	
  	char DSSP_Seq[1000],SS[1000],PATH[1000];		  
	int matrix[1001][1001];
	float rmsd;
	float cutoff;
	FILE *in;

	if(argc<2)
	{
		printf("USAGE: %s LIB_FILE\n",argv[0]);
		return 0;
	}

  	/*** FILE HANDLING: ***/
	in = fopen(argv[1],"r");
	
	for(k=0;fscanf(in,"%s\t%c\t%d\t%d\t",Header,&Chain,&start1,&length)!=EOF ;k++)
	{
		
		//fscanf(in,"%s\t%c\t%d\t%d\t",Header,&Chain,&start1,&length)

		if ( (Atmrec.frag[0]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}

		fscanf(in,"%s\t%c",Atmrec.frag[0]->ALN_Seq,&Atmrec.frag[0]->aux);
		fscanf(in,"\t%d\t%d\t%d\t%d\t%lf\t%d\t%lf\n",&Atmrec.frag[0]->match_score,&Atmrec.frag[0]->seq_score,&LLL.length,&Atmrec.frag[0]->start2,&Atmrec.frag[0]->resolution,&Atmrec.frag[0]->ss_score,&Atmrec.frag[0]->rmsd);
		strcpy(Atmrec.frag[0]->Header,Header);
		Atmrec.frag[0]-> Chain = Chain;
		Atmrec.frag[0]-> start_res	= start1;
		Atmrec.frag[0]-> end_res	= length;
		Atmrec.frag[0]-> length	= LLL.length;
		sprintf(Atmrec.frag[0]->fname,"%s/%c%c%c%c.pdb",getenv("PDB"),tolower(Header[0]),tolower(Header[1]),tolower(Header[2]),tolower(Header[3]));
		if(Atmrec.frag[0]->rmsd<1.5)
		{
			extract_frags(Atmrec.frag[0]->fname,LLL.length,Atmrec.frag[0]->start_res,Chain,Atmrec.frag[0]->ALN_Seq);
			sprintf(AUX2,"mv ./temp.pdb ./1AILfrag%d.pdb",numfrags);
			printf("hide lines, 1AILfrag%d\nshow cartoon, 1AILfrag%d\n",numfrags,numfrags);
			printf("align 1AILfrag%d, 1AIL and resi %d-%d\n",numfrags,Atmrec.frag[0]->start2,Atmrec.frag[0]->start2+LLL.length);
			numfrags++;
			system(AUX2);
			
		}
	}

	fclose(in);
	return 0;  
}


