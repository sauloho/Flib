#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>     

#define MATCH 2
#define MISMATCH -2
#define score(A,B) (((A)==(B))?(MATCH):(MISMATCH))

int Blossum[24][24];          /* The Env. Matrix (CHANGE NAME WHEN NOT BLOSSUM).    */
char Dict[24]={'A','R','N','D','C','Q','E','G','H','I','L',
               'K','M','F','P','S','T','W','Y','V','B','Z',
               'X','*'};      /* A one-letter code residue dictionary.              */

int find(char A)
{
  int i;
  for(i=0;i<24;i++)
    if(Dict[i]==A) 
      return i; 
  return 23;  
}

/* Sequence Alignment scoring function */
int score3(char A, char B) 
{
  return Blossum[find(A)][find(B)];  
}


int score2(char A, char B) 
{
	if( A=='I' || A =='T' || A=='S' )
		A=' ';
	if( B =='I' || B=='T' || B=='S' )
		B=' ';

	if( A == 'G' && B == ' ')
		return 1;
	if(A=='G' && B=='H')
		return 1;
	if(B=='G' && A==' ')
		return 1;
	if(B=='G' && A=='H')
		return 1;
	if(A=='B' && B==' ')
		return 1;
	if(A=='B' && B=='E')
		return 1;
	if(B=='B' && A==' ')
		return 1;
	if(B=='B' && A=='E')
		return 1;

	return score(A,B);
}


int main(int argc, char *argv[])
{
	int length,start1,match_score,seq_score,start2,ss_score,total,old=-1;
	int i,loop,beta,helix,j;
	char Header[82],Chain,type,ALN_Seq[200];
	char Fasta_Seq[1000];
	char AUX[82],Seq2[10],SS1[10],SS2[10];
	char PATH[1000];
	float rmsd;
	double resolution;
	FILE *library,*blossum_file;

	if(!argv[1])
	{
		printf("USAGE: %s PDB_ID.rmsd_lib PDB_ID\n",argv[0]);
		return 0;
	}

  	/*** FILE HANDLING: ***/
	library = fopen(argv[1],"r");
        sprintf(PATH,"%s/data/blossum62.txt",getenv("FLIB"));
	blossum_file = fopen(PATH,"r");

	/* READ THE ENV. MATRIX FOR THE ALIGNMENT */
	for(i=0;i<24;i++)
	    	for(j=0;j<24;j++)
				fscanf(blossum_file,"%d",&Blossum[i][j]); 

	for(total=0;fscanf(library,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%s\t%s\t%s\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,Seq2,SS1,SS2)!=EOF;total++)
	{
		loop=0; beta=0; helix=0;
		for(i=0;i<9;i++)
		{
			switch(SS2[i])
			{
				case 'H':
					helix++;
					break; 
				case 'E':
					beta++;
					break;
				default:
					loop++;
					break; 
			}
			ss_score+=score2(SS1[i],SS2[i]);
			seq_score+=score3(Seq2[i],ALN_Seq[i]);	
		}
		if (helix>=5)
			type='H';
		else if(beta>=5)
			type='B';
		else if(loop>=5)
			type='L';
		else
			type='O';
		printf("%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t0.00\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score);
	}

	fclose(library);
	return 0;
}
