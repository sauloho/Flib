 /*******************************************************************/
/* 								   */		
/*     C. Deane - Superposition Function - 27/10/04		   */
/* 								   */
/*	Take in lists of fragments perform pariwise superpositon.  */
/*	then cluster on the basis of these distances into families.*/
/* 								   */
/*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "slink.h"

#define MaxFrag 10000
#define TRUE 1
#define FALSE 0

#define min(a,b) (((a)<(b))?(a):(b))

BITS Atmrec;
int num_printed;
int Vector[1001];

/* Collect the coordinates from the pdb file starting from the correct residue. */
int coordcol( char bum[FILEL],int u, int length, int start_res,char ChainQ,char* Seq);
/*perform the superposition*/
float super4(FRAGMENTS *One, FRAGMENTS *Two, int nres);

int main (int argc,char *argv[])
{
	/* this is a list of file names*/
	int k,total,status,temp_nres,temp_natom;
	int old_res=-1;
	int frag;
	float rmsd;
	float cutoff;
	FILE *in;

	if(argc<3)
	{
		printf("USAGE: %s PDB_ID PATH_TO_PDB cutoff\n",argv[0]);
		return 0;
	}

  	/*** FILE HANDLING: ***/
	in = fopen(argv[1],"r");
	cutoff=atof(argv[3]);
	frag=atoi(argv[4]);

	if ( (Atmrec.frag[0]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}
	if ( (Atmrec.frag[1]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}
	
	for(k=0;fscanf(in,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%lf\n",Atmrec.frag[1]->Header,&Atmrec.frag[1]->Chain,&Atmrec.frag[1]->start_res,&Atmrec.frag[1]->end_res,Atmrec.frag[1]->ALN_Seq,&Atmrec.frag[1]->aux,&Atmrec.frag[1]->match_score,&Atmrec.frag[1]->seq_score,&		Atmrec.frag[1]->length,&Atmrec.frag[1]->start2,&Atmrec.frag[1]->resolution,&Atmrec.frag[1]->ss_score,&Atmrec.frag[1]->rmsd)!=EOF ;k++)
	{
		if(Atmrec.frag[1]->start2 != old_res)	/* This is a new residue! */
		{
			k=0;
			total=0;
			old_res = Atmrec.frag[1]->start2;
		}
		if(k==frag)
		{
			strcpy(Atmrec.frag[0]->Header,Atmrec.frag[1]->Header);
			Atmrec.frag[0]->Chain = Atmrec.frag[1]->Chain;
			Atmrec.frag[0]->start_res = Atmrec.frag[1]->start_res;
			Atmrec.frag[0]->end_res = Atmrec.frag[1]->end_res;
			strcpy(Atmrec.frag[0]->ALN_Seq,Atmrec.frag[1]->ALN_Seq);
			Atmrec.frag[0]->aux = Atmrec.frag[1]->aux;
			Atmrec.frag[0]->match_score = Atmrec.frag[1]->match_score;
			Atmrec.frag[0]->seq_score = Atmrec.frag[1]->seq_score;
			Atmrec.frag[0]->length = Atmrec.frag[1]->length; 
			Atmrec.frag[0]->start2 = Atmrec.frag[1]->start2;
			Atmrec.frag[0]->resolution = Atmrec.frag[1]->resolution;
			Atmrec.frag[0]->ss_score = Atmrec.frag[1]->ss_score;
			Atmrec.frag[0]->rmsd = Atmrec.frag[1]->rmsd;
			
			sprintf(Atmrec.frag[0]->fname,"%s/%c%c%c%c.pdb",argv[2],tolower(Atmrec.frag[0]->Header[0]),tolower(Atmrec.frag[0]->Header[1]),tolower(Atmrec.frag[0]->Header[2]),tolower(Atmrec.frag[0]->Header[3]));
			status = coordcol(Atmrec.frag[0]->fname,0,Atmrec.frag[0]->length,Atmrec.frag[0]->start_res,Atmrec.frag[0]->Chain,Atmrec.frag[0]->ALN_Seq);
			
			if(!status) k--;
		}
		if (k >= 9 && total < 10 )
		{
			sprintf(Atmrec.frag[1]->fname,"%s/%c%c%c%c.pdb",argv[2],tolower(Atmrec.frag[1]->Header[0]),tolower(Atmrec.frag[1]->Header[1]),tolower(Atmrec.frag[1]->Header[2]),tolower(Atmrec.frag[1]->Header[3]));
			status = coordcol(Atmrec.frag[1]->fname,1,Atmrec.frag[1]->length,Atmrec.frag[1]->start_res,Atmrec.frag[1]->Chain,Atmrec.frag[1]->ALN_Seq);
			if(status==1)
			{
				fprintf(stderr,"Could not load fragment %d for position %d!\n",k,Atmrec.frag[1]->start2);
				continue;
			}

			status = 0;
			if(Atmrec.frag[0]->numres != Atmrec.frag[1]->numres && Atmrec.frag[0]->numres && Atmrec.frag[1]->numres )
			{	
				if(Atmrec.frag[0]->numres > Atmrec.frag[1]->numres)
				{
					status=1;
					temp_nres = Atmrec.frag[0]->numres;
					temp_natom = Atmrec.frag[0]->natom;
					Atmrec.frag[0]->numres = Atmrec.frag[1]->numres;
					Atmrec.frag[0]->natom = Atmrec.frag[1]->natom;
				}
				else
				{
					status=2;
					temp_nres = Atmrec.frag[1]->numres;
					temp_natom = Atmrec.frag[1]->natom;
					Atmrec.frag[1]->numres = Atmrec.frag[0]->numres;
					Atmrec.frag[1]->natom = Atmrec.frag[0]->natom;
				}
			}
			rmsd = super4(Atmrec.frag[0],Atmrec.frag[1],Atmrec.frag[1]->numres);
	
			if(rmsd < cutoff && rmsd > -0.1)
			{	
				total++;
				printf("%s\t%c\t%3d\t%3d\t",Atmrec.frag[1]->Header,Atmrec.frag[1]->Chain,Atmrec.frag[1]->start_res,Atmrec.frag[1]->end_res);
				printf("%s\t%c",Atmrec.frag[1]->ALN_Seq,Atmrec.frag[1]->aux);
				printf("\t%d\t%d\t%d\t%3d\t%.3lf\t%d\t%.3lf\n",Atmrec.frag[1]->match_score,Vector[1],Atmrec.frag[1]->length,Atmrec.frag[1]->start2,Atmrec.frag[1]->resolution,Atmrec.frag[1]->ss_score,Atmrec.frag[1]->rmsd);
			}

			if(status==1)
			{
				Atmrec.frag[0]->numres = temp_nres;
				Atmrec.frag[0]->natom = temp_natom;
			}
			if(status==2)
			{
				Atmrec.frag[1]->numres = temp_nres;
				Atmrec.frag[1]->natom = temp_natom;
			}	
		
		}
	}
	fclose(in);
	return 0;  
}

