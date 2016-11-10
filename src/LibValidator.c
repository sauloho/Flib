#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<ctype.h>
#include"slink.h"

#define MAXLEN 10000
#define MAX_FRAG 3000

/* Collect the coordinates from the pdb file starting from the correct residue. */
int coordcol( char bum[FILEL],int u, int length, int start_res,char ChainQ,char* Seq);
/* Performs the fragment superposition, returns the RMSD */
float super4(FRAGMENTS *One, FRAGMENTS *Two, int nres);
BITS Atmrec;
FRAGMENTS *all;

/*typedef struct Fragment
{
        char line[250]; 
        struct Fragment *next;
}fragment;*/

//safdgasfgas
int main(int argc, char* argv[])
{
        int length,seq_score,ramach_score,ss_score,tor_score,top=0;
        int i,j,k,start1,start2,start_res,m; 							/* Counters           				       */
        char AUX[600],SS[2],Res[2];       							/* Multi-purpose Auxiliary Strings         */
        char c;				                                        /* Auxiliary Char                          */
        char DB_Seq[MAXLEN],Fasta_Seq[MAXLEN];              /* The DB protein fasta sequence                         */
        char Chain;				                  /* The Protein chain.                                        */
        char Header[82];        			      /* Information extracted from the header on the Database */
        float rmsd;
        double resolution;                        /* The resolution of the crystal structure on the DB.    */
        double new_score;                         /* The predicted torsion angle score                                     */

        FILE *library, *input_fasta;

        if(argc < 1 )
        {
                printf("USAGE: %s PDB LIBRARY.lib Chain\n",argv[0]);
                return 0;
        }

        /*** FILE HANDLING: ***/

        strcpy(AUX,argv[1]);
        input_fasta = fopen(strcat(AUX,".fasta.txt"),"r");
        if (input_fasta == NULL) {printf("Fasta file not found: %s\n", AUX); return 0;}

        library = fopen(argv[2],"r");
        if(library==NULL) { fprintf(stderr,"Library file: %s not found!\n",argv[2]);    return 0;       }

		
	if ( (Atmrec.frag[0]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}
        if ( (Atmrec.frag[1]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}

        /***** READ QUERY FASTA SEQUENCE *****/

        /* Remove Header from FASTA file */
        for(c = fgetc(input_fasta); c!='\n' ; c=fgetc(input_fasta));
        /* Read the sequence and store it in the Fasta_Seq string */
        for(fscanf(input_fasta,"%s",Fasta_Seq); fscanf(input_fasta,"%s",AUX)!=EOF ; strcat(Fasta_Seq,AUX) );
        /* Compute the length of the sequence and store it in m */
        m = strlen(Fasta_Seq);
        /***** END OF READ QUERY FASTA SEQUENCE *****/



	/*** FRAGMENT VALIDATION ****/
	while(  fscanf(library,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%lf",Header,&Chain,&start1,&k,DB_Seq,&c,&ramach_score,&seq_score,&length,&start2,&resolution,&ss_score,&new_score) != EOF )
	{

        	top=0;
	        Atmrec.frag[0]->start_res       = start1;
         	Atmrec.frag[1]->start_res       = start2;

	        sprintf(Atmrec.frag[0]->fname,"%s/%c%c%c%c.pdb",getenv("PDB"),tolower(Header[0]),tolower(Header[1]),tolower(Header[2]),tolower(Header[3]));
	        top = coordcol(Atmrec.frag[top]->fname,top,length,start1,Chain,DB_Seq);
		if(top<1)
		{
                        fprintf(stderr,"Could not open the fragment starting at position %d for %s_%c.\n",start1,Header,Chain);
			continue;
		}
	        sprintf(Atmrec.frag[1]->fname,"./%s.pdb",argv[1]);                                   
		for (i=0;i<length;i++)
			AUX[i]=Fasta_Seq[start2+i];
		AUX[i]='\0'; 
	        top = coordcol(Atmrec.frag[top]->fname,top,length,start2,argv[3][0],AUX);
  		/* If, by any reason, was not capable of finding the fragment: */
	        if(top<2)
		{
                        fprintf(stderr,"Unresolved fragment at Xstal structure at position %d (length: %d).\n",start2,length);
                        printf("%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t-1.0\n",Header,Chain,start1,k,DB_Seq,c,ramach_score,seq_score,length,start2,resolution,ss_score,new_score);

		}
         	else
		{
			/* Sanity check to see if they have the same number of atoms */
 			if( Atmrec.frag[0]->natom == Atmrec.frag[1]->natom  &&  Atmrec.frag[1]->natom > 0  &&  Atmrec.frag[0]->natom > 0 && Atmrec.frag[0]->numres > 0 && Atmrec.frag[1]->numres > 0 && Atmrec.frag[0]->numres == Atmrec.frag[1]->numres )
	        	{
				rmsd =  super4( Atmrec.frag[0], Atmrec.frag[1],length);
		                if(rmsd > 10.0 || rmsd < 0.0)
				{
					fprintf(stderr,"Could not align the fragment starting at position %d for %s_%c.\n",start1,Header,Chain);
					continue;
				}
				else
		                    printf("%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\t%.2f\n",Header,Chain,start1,k,DB_Seq,c,ramach_score,seq_score,length,start2,resolution,ss_score,new_score,rmsd);
			}
        	}
	}


	fclose(library);
	return 0;
}
