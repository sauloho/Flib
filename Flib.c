#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<ctype.h>
#include"slink.h"

#define RANDOM 1
#define EXHAUSTIVE 1
#define VALIDATE 0
#define MaxFrag 10000
#define TRUE 1
#define FALSE 0
#define CUTOFF 1.5
#define MATCH 2
#define MISMATCH -2
#define GAP -100
#define VERBOSE 1
#define MAXLEN 10000
#define MAX_FRAG 3000

#define max(A,B) (((A)>(B))?(A):(B))
#define min(A,B) (((A)<(B))?(A):(B))
#define score(A,B) (((A)==(B))?(MATCH):(MISMATCH))

/* Collect the coordinates from the pdb file starting from the correct residue. */
int coordcol( char bum[FILEL],int u, int length, int start_res,char ChainQ,char* Seq);
int coordcol2( FRAGMENTS * Piece,char PDB[82],char Chain, char* Seq,int *start );

/* Performs the fragment superposition, returns the RMSD */
float super4(FRAGMENTS *One, FRAGMENTS *Two, int nres);
float Calc_dih(RESIDUE * OneP,RESIDUE * One,RESIDUE * OneA, float * temp);
BITS Atmrec;
FRAGMENTS *all;

typedef struct confidence
{
	float coil;
	float helix;
	float beta;
}CONF;

typedef struct Fragment 
{
	char line[250];
	struct Fragment *next;
}fragment;

int   A[600][MAXLEN];			  /* The Alignment Matrix 								   */
int   A_len[600][MAXLEN];		  /* The Alignment Length Matrix						   */
int   A_SS[600][MAXLEN];		  /* The Alignment SS score Matrix						   */
int   A_ramach[600][MAXLEN];	  /* The Alignment Ramach score Matrix	 				   */
float A_SSconf[600][MAXLEN];	  /* The Alignment SS Conf score Matrix					   */
fragment *Frag_lib[1000][100]; 	  /* A linked list for each query position (max=1000). Each sequence score is a key for the list (max=100). */
fragment *Frag_lib_rnd[1000][100]; 	  /* A linked list for each query position (max=1000). Each sequence score is a key for the list (max=100). */

/* Given a one-letter code for an aa, returns the index */
/* corresponding to that aa on the dictionary.          */
int find(char A,char *Dict)
{
  int i;
  for(i=0;i<strlen(Dict);i++)
    if(Dict[i]==A) 
      return i; 
  return strlen(Dict)+1;  
}

/* SS Alignment scoring function */
int score_ss(char A, char B) 
{
	if( A=='C' )
		A=' ';
	if( B =='I' || B=='T' || B=='S' )
		B=' ';

	if(B=='G' && A==' ')
		return 1;
	if(B=='G' && A=='H')
		return 1;
	if(B=='B' && A==' ')
		return 1;
	if(B=='B' && A=='E')
		return 1;

	return score(A,B);
}

float score_ssconf(CONF confidence, char B) 
{

	if( B!='H' && B !='E' )
		return confidence.coil;
	if( B=='H' )
		return confidence.helix;
	if( B=='E' )
		return confidence.beta;
	return 0.0;
}

int main(int argc,char* argv[])
{
	short int LOGODS[7][20][20];  /* FREAD Environment Tables                              */
	short int Blossum[24][24];    /* BLOSSUM62 Env. Table                                  */
	int length,seq_score,ramach_score,ss_score,top=0;
	int i,j,k,k2,total,start1,start2,index1,index2,start_res; /* Counters                  */
	int loop,helix,beta;		  /* Counters											   */
	int *Total;			 		  /* Total # of frags for each query position.	  		   */
	int *Total_rnd;		 		  /* Total # of frags for each query position.	  		   */
	int *Worst;					  /* Score of worst frag for each query position.		   */	
	int *Worst_rnd;				  /* Score of worst frag for each query position.		   */	
	int m;			 			  /* The length of the Query's fasta sequence              */
	int n;			  		      /* The length of the DB protein fasta sequence           */
	float Phi[MAXLEN], Psi[MAXLEN];
	char AUX[600],SS[2],Res[2];	  /* Multi-purpose Auxiliary Strings 		               */
	char c;						  /* Auxiliary Char 						               */
	char DB_Seq[MAXLEN];		  /* The DB protein fasta sequence                         */
  	char DB_SS[MAXLEN];		      /* The DB protein's secondary structure        		   */
	int DB_Ramach[MAXLEN];		  /* The DB protein's Ramachandran grouping.		       */	
	char Fasta_Seq[MAXLEN];		  /* The Query's fasta sequence                            */
	char Fasta_SS[MAXLEN];		  /* The Query's predicted secondary structure             */
	char Chain;                   /* The Protein chain.                            		   */
	char Header[82];              /* Information extracted from the header on the Database */
	char Query[82];
	char Dict2[21]={'G','A','V','L','M','I','F','Y','W','S','T','C','P','N','Q','K','R','H','D','E','\0'};              /* A one-letter code residue dictionary. */	
	char PATH[1000];
	int path_len;
	float rmsd;
	float probability;
	float Angles[MAXLEN][4];
	double resolution;			  /* The resolution of the crystal structure on the DB.	   */	
	double new_score;			  /* The predicted torsion angle score					   */	

	CONF Fasta_Conf[MAXLEN];      /* The confidence in the predicted secondary structure   */
	fragment *new_frag;
	FILE *input_fasta,*input_ss,*input_pdb,*input_phipsi,*blossum_file,*logods;
	srand (time(NULL));

	if(argc != 2  )
	{
		printf("Usage: %s PDB_CODE \n",argv[0]);
		return 0;
	}

	strcpy(Query,argv[1]);


  	/*****      FILE HANDLING   *****/
	strcpy(AUX,argv[1]);
	input_fasta = fopen(strcat(AUX,".fasta.txt"),"r");
	if (input_fasta == NULL) {printf("Fasta file not found: %s\n", AUX); return 0;}	
	
	strcpy(AUX,argv[1]);
	input_ss = fopen(strcat(AUX,".fasta.ss"),"r");
	if (input_ss == NULL) {printf("Predicted Secondary Structure file not found: %s\n", AUX); return 0;}

	sprintf(PATH,"%s/parsedPDB_new.txt",getenv("FLIB"));
	input_pdb = fopen(PATH,"r");
	if (input_pdb == NULL) {printf("Protein database file is missing: parsedPDB_new.txt\n"); return 0;}

	strcpy(AUX,argv[1]);
	input_phipsi = fopen(strcat(AUX,".spXout"),"r");
	if (input_phipsi == NULL) {printf("SPINE-X input file (predicted torsion angles) not found!\n"); return 0;}

        sprintf(PATH,"%s/blossum62.txt",getenv("FLIB"));
	blossum_file = fopen(PATH,"r");
	if (blossum_file == NULL) {printf("BLOSSUM file not found!\n"); return 0;}


	/***** END OF FILE HANDLING *****/

	/***** READ QUERY FASTA SEQUENCE *****/

	/* Remove Header from FASTA file */
	for(c = fgetc(input_fasta); c!='\n' ; c=fgetc(input_fasta)); 
	/* Read the sequence and store it in the Fasta_Seq string */
	for(fscanf(input_fasta,"%s",Fasta_Seq); fscanf(input_fasta,"%s",AUX)!=EOF ; strcat(Fasta_Seq,AUX) );
	/* Compute the length of the sequence and store it in m */
	m = strlen(Fasta_Seq);
	/***** END OF READ QUERY FASTA SEQUENCE *****/

	/***** READ QUERY'S PREDICTED SECONDARY STRUCTURE SEQUENCE *****/

	/* READ INPUT SECONDARY STRUCTURE SEQUENCE */
	while(fscanf(input_ss,"%d",&i)!=EOF && i<= m  )
		fscanf(input_ss,"%s %c %f %f %f",AUX,&Fasta_SS[i-1],&Fasta_Conf[i-1].coil,&Fasta_Conf[i-1].helix,&Fasta_Conf[i-1].beta);
	if(i!=m) { printf("ERROR: Fasta sequence and pred. secondary structure have different lengths!\n"); return 0; }
	Fasta_SS[m]='\0'; 
	/***** END OF READ QUERY'S PREDICTED SECONDARY STRUCTURE SEQUENCE *****/
	
	/***** READ QUERY'S PREDICTED TORSION ANGLES *****/
	
	/* Remove Header from SPINE-X output file */
	for(c=fgetc(input_phipsi);c!='\n';c=fgetc(input_phipsi));

	for(i=0;fscanf(input_phipsi,"%d %s %s %f %f",&k,Res,SS,&Angles[i][0],&Angles[i][1])!=EOF;i++)
	{
		Angles[i][2]=1.0; Angles[i][3]=1.0;
		for(c=fgetc(input_phipsi);c!='\n';c=fgetc(input_phipsi));
	}

	/***** END OF READ QUERY'S PREDICTED TORSION ANGLES *****/
	
	/***** READ BLOSSUM MATRIX FOR THE ALIGNMENT *****/
	for(i=0;i<24;i++)
		for(j=0;j<24;j++)
			fscanf(blossum_file,"%d",&Blossum[i][j]);  
	/***** END OF READ BLOSSUM MATRIX FOR THE ALIGNMENT *****/
	
	/***** READ FREAD ENV. TABLES *****/
	sprintf(AUX,"%s/logods0.txt",getenv("FLIB"));
	path_len = strlen(AUX);
	for(;AUX[path_len-5]<48+7;AUX[path_len-5]++)
	{	
		logods=fopen(AUX,"r");
		if (logods == NULL) {printf("%s file for FREAD environment matrix not found\n",AUX); return 0;}

		for(i=0;i<20;i++)
			for(j=0;j<20;j++)
				fscanf(logods,"%hd",&LOGODS[(int)AUX[path_len-5]-48][i][j]);
		fclose(logods);
	}
	/***** END OF READ FREAD ENV. TABLES *****/
	
	/***** INITIALIZE THE FRAGMENT LIBRARY LINKED LIST *****/
	Total=malloc(sizeof(int)*m);
	Total_rnd=malloc(sizeof(int)*m);
	Worst=malloc(sizeof(int)*m);
	Worst_rnd=malloc(sizeof(int)*m);
	for(i=0;i<m;i++)
	{
		Total[i]=0;
		Total_rnd[i]=0;
		Worst[i]=0;
		Worst_rnd[i]=-1;
		for(j=0;j<100;j++)
		{
			Frag_lib[i][j]=NULL;
			Frag_lib_rnd[i][j]=NULL;
		}
	}

	/***** END OF INITIALIZE THE FRAGMENT LIBRARY LINKED LIST *****/

	/***** BEGINNING OF FRAGMENT EXTRACTION *****/

	probability = (3.0*(float)(m)*10000)/(120916.0/4); /* Define a probability for random fragment extraction */
	/* Iterate over every sequence in the input database */
	for(total=0; fscanf(input_pdb,"%s %s %c %lf",AUX,Header,&Chain,&resolution) != EOF  ;total++)
	{
		/* Read DB protein fasta sequence */
		fscanf(input_pdb,"%s",DB_Seq);			/* Extract the protein sequence from the database. */	
		fgetc(input_pdb); 					    /* Get the line break char (chomp) 				   */
		n=strlen(DB_Seq);						/* Get the length of the DB protein sequence 	   */

		//printf("%s %d\n",DB_Seq,n);
		
		/* Read DB protein secondary structure  */
		for(j=0 , DB_SS[j]= fgetc(input_pdb); DB_SS[j]!=(int)'\n' && DB_SS[j] != EOF ; j++, DB_SS[j]= fgetc(input_pdb));
		DB_SS[j]='\0';
		if(j != n)
		{
			for(c=fgetc(input_pdb);c!='>';c=fgetc(input_pdb));
			ungetc(c,input_pdb);
			fprintf(stderr,"Sequence and secondary structure uncompatible lengths for protein: %s\t Seq: %d\tSS: %d\n",Header,n,j);
			continue;
		}
		
		/* Read DB protein Ramachandran Groups */
		for(j=0 , DB_Ramach[j]= (int)fgetc(input_pdb); DB_Ramach[j]!=(int)'\n' && DB_Ramach[j] != EOF ; j++, DB_Ramach[j]= (int)fgetc(input_pdb))
				DB_Ramach[j]-=48;	
		if(j!= n)
		{
			for(c=fgetc(input_pdb);c!='>';c=fgetc(input_pdb));
			ungetc(c,input_pdb);
			fprintf(stderr,"Sequence and Ramachandran Grouping uncompatible lengths for protein: %s\t Seq: %d\tRamach: %d\n",Header,n,j);
			continue;
		}

		/***** RANDOM EXTRACTION *****/
		#if RANDOM
		for(j=probability;j >= 0 && n >= 15; j--)
		{	
			if(!j && (float)(rand()%10000)/10000.0 > probability-(int)probability)
		    		continue;

			length = rand() % 7 + 6; /* Randomize fragment length between 6 and 12 residues */
			start1  = rand() % (m-length + 1); /* Randomize fragment starting position in query between 0 and (length of protein - length of fragment) */
			start2  = rand() % (n-length + 1); /* Randomize fragment starting position in pdb structure between 0 and (length of protein - length of fragment) */
			loop=0; helix=0; beta=0; ss_score=0; seq_score=0; ramach_score=0;

			/* Compute sequence and secondary structure score for the fragment */
			for (k=start1,k2=start2; k-start1 < length; k++,k2++)
			{
				index1=find(Fasta_Seq[k],Dict2);
				index2=find(DB_Seq[k2],Dict2);
				if((DB_Ramach[start2] >= 7 || DB_Ramach[start2] < 0 ) || (index1 < 0 || index1 >= 20 || index2 < 0 || index2 >=20 ))
				{
					ramach_score=-50;
					ss_score=-50;
					break;
				}
				else
					ramach_score += LOGODS[DB_Ramach[start2]][index1][index2];
		       		ss_score += score_ss(Fasta_SS[k],DB_SS[k2]);
				/* Find predominant SS on the fragment */
				switch(DB_SS[k2])
				{
					case ' ':
						loop++;
						break;
					case 'H':
						helix++;
						break; 
					case 'E':
						beta++;
						break;
					default:
						break; 
				}
			}
			
			if( ( ramach_score > 6 || ( helix < length/2+1 && beta < length/2 + 1 && ramach_score > 0 )) && (Total_rnd[start1] < MAX_FRAG || ramach_score > Worst_rnd[start1] ) )
			{	
				if(!( helix < length/2+1 && beta < length/2 + 1 ) && ss_score < 2 )
					continue; 

				new_frag=malloc(sizeof(fragment));
				if(new_frag == NULL)
				{
					fprintf(stderr,"Malloc error!\n");
					break;
				}
				if(ramach_score > 99) ramach_score = 99; /* This should not happen very often! */

				new_frag->next = Frag_lib_rnd[start1][ramach_score];

				sprintf(new_frag->line,"%s\t%c\t%3d\t%3d\t",Header,Chain,start2,start2+length-1);

				for (k=0; k < length ; k++)
					AUX[k]=DB_Seq[start2+k];
				AUX[k]='\0';
				strcat(new_frag->line,AUX);

				if(helix>length/2+1)
					strcat(new_frag->line,"\tH");
				else { 	if(beta>length/2+1)
					strcat(new_frag->line,"\tB");
				else { 	if(loop>length/2+1)
					strcat(new_frag->line,"\tL");
				else 
					strcat(new_frag->line,"\tO"); } }

				sprintf(AUX,"\t%d\t%d\t%d\t%d\t%.2lf\t%d\n",ramach_score,0,length,start1,resolution,ss_score);
				strcat(new_frag->line,AUX);
				Frag_lib_rnd[start1][ramach_score] = new_frag;

				if(Total_rnd[start1]==0) /* First fragment for that position! */
						Worst_rnd[start1]=ramach_score; /* Update the score of the Worst Fragment for that position */	

				if(Total_rnd[start1] < MAX_FRAG)
				{
					Total_rnd[start1]++;	/* Update the number of fragments for that position. */
					if(ramach_score < Worst_rnd[start1]) Worst_rnd[start1]=ramach_score; /* Update the score of the Worst Fragment for that position */
				}
				else 
				{
					/* Remove the worst one */
					if(Frag_lib_rnd[start1][Worst_rnd[start1]] != NULL)
					{
						new_frag = Frag_lib_rnd[start1][Worst_rnd[start1]];
						Frag_lib_rnd[start1][Worst_rnd[start1]] = new_frag->next;	
						free(new_frag);
					}
					if( Frag_lib_rnd[start1][Worst_rnd[start1]] == NULL )
						for( ; Worst_rnd[start1] < 100; Worst_rnd[start1]++)
							if(Frag_lib_rnd[start1][Worst[start1]]!=NULL)
								break;
				}	

			}
		}
		#endif

		/***** EXHAUSTIVE EXTRACTION *****/
		#if EXHAUSTIVE	
		/***** PERFORM THE SEQUENCE ALIGNMENT *****/		

		/* Initialize the Scoring Matrix */
		for(i=0;i<=m;i++) {
			A[i][0]=0; 		A_len[i][0]=0; A_SS[i][0]=0; A_SSconf[i][0]=0; A_ramach[i][0]=0;	}
		for(j=0;j<=n;j++) {	
			A[0][j]=0;		A_len[0][j]=0; A_SS[0][j]=0; A_SSconf[0][j]=0; A_ramach[0][j]=0;	}	

		/* Populate the scoring matrix using sequence score! */

		for(i=1;i<=m;i++)
		{
			for(j=1;j<=n;j++)
    		{
				
   				if(DB_Ramach[j-1] >= 0 && DB_Ramach[j-1] < 7 )
				{
					index1=find(Fasta_Seq[i-1],Dict2);
					index2=find(DB_Seq[j-1],Dict2);
					if(index1 < 0 || index1 >= 20 || index2 < 0 || index2 >=20 )
					 	A_ramach[i][j]=0;
					else
						A_ramach[i][j] = max( A_ramach[i-1][j-1]+LOGODS[DB_Ramach[j-1]][index1][index2] , 0 );
				}
				else
					A_ramach[i][j] = 0;
				
				if( !A_ramach[i][j] )
				{
 					 A[i][j]       = 0;
					 A_len[i][j]   = 0;
					 A_SS[i][j]    = 0;
					 A_SSconf[i][j]= 0.0;

				}
				else
				{
		 		//	 A[i][j]       = A[i-1][j-1]       + Blossum[find(Fasta_Seq[i-1],Dict1)][find(DB_Seq[j-1],Dict1)];
					 A[i][j]	   = 0;
					 A_len[i][j]   = A_len[i-1][j-1]   + 1;
					 A_SS[i][j]    = A_SS[i-1][j-1]    + score_ss(Fasta_SS[i-1],DB_SS[j-1]);
					 A_SSconf[i][j]= A_SSconf[i-1][j-1]+ score_ssconf(Fasta_Conf[i-1],DB_SS[j-1]);
				}
			} 
		}


		/* "Traceback" step of sequence alignment */   
		for(i=m;i>0;i--)
			for(j=n;j>0;j--)
			{

				/* If this fragment has a better score than the worst fragment we extracted so far 
				   for that position and if the lenght of the fragment is greater than 5, we extract it!     */
				if(A_len[i][j] > 5)
				{
					start1 = i - A_len[i][j];
					start2 = j - A_len[i][j];
				}
				else
					continue;

				if( start1 < 0 || start1 >= m || start2 < 0 || start2 >=n ) continue;				
							
				if( (Total[start1] < MAX_FRAG || A_ramach[i][j] > Worst[start1] )  && A_len[i][j] > 5 /*&& A_SS[i][j] > 0*/ && A_len[i][j] < 20 /*A_ramach[i][j] > 50*/)
				{	
					loop=0; helix=0; beta=0;
					/* Find predominant SS on the fragment */
					for (k=0; k < A_len[i][j] && k + start2 < n; k++)
					{
						AUX[k]=DB_Seq[start2+k];
						switch(DB_SS[start2+k])
						{
								case ' ':
									loop++;
									break;
								case 'H':
									helix++;
									break; 
								case 'E':
									beta++;
									break;
								default:
									break; 
						}
					}
					AUX[k]='\0';

					if( !( helix < A_len[i][j]/2+1 && beta < A_len[i][j]/2 + 1 ) && A_SS[i][j] < 2 && A_ramach[i][j] < 70  )
						continue;


					if(A_ramach[i][j] > 99) A_ramach[i][j] = 99; /* This should not happen very often! */

					/* Allocate the fragment */
					new_frag = malloc(sizeof(fragment));
					new_frag->next = Frag_lib[start1][A_ramach[i][j]];
					sprintf(new_frag->line,"%s\t%c\t%3d\t%3d\t",Header,Chain,start2,start2+A_len[i][j]-1);


					strcat(new_frag->line,AUX);
					if(helix>=A_len[i][j]/2+1)
						strcat(new_frag->line,"\tH");
					else { 	if(beta>=A_len[i][j]/2+1)
						strcat(new_frag->line,"\tB");
					else { 	if(loop>=A_len[i][j]/2+1)
						strcat(new_frag->line,"\tL");
					else 
						strcat(new_frag->line,"\tO"); } }

					sprintf(AUX,"\t%d\t%d\t%d\t%d\t%.2lf\t%d\n",A_ramach[i][j],A[i][j],A_len[i][j],start1,resolution, A_SS[i][j]);
					strcat(new_frag->line,AUX);

					Frag_lib[start1][A_ramach[i][j]] = new_frag;
					/* If the number of fragments for that position is less than MAX_FRAG... */

					if(Total[start1]==0) /* First fragment for that position! */
							Worst[start1]=A_ramach[i][j]; /* Update the score of the Worst Fragment for that position */	

					if(Total[start1]<MAX_FRAG)
					{
						Total[start1]++;	/* Update the number of fragments for that position. */
						if(A_ramach[i][j] < Worst[start1]) Worst[start1]=A_ramach[i][j]; /* Update the score of the Worst Fragment for that position */
					}
					else 
					{
							/* Remove the worst one */
							if(Frag_lib[start1][Worst[start1]] != NULL)
							{
								new_frag = Frag_lib[start1][Worst[start1]];
								Frag_lib[start1][Worst[start1]] = Frag_lib[start1][Worst[start1]]->next;	
								free(new_frag);
							}
							if( Frag_lib[start1][Worst[start1]] == NULL )
								for( ; Worst[start1] < 100; Worst[start1]++)
									if(Frag_lib[start1][Worst[start1]]!=NULL)
										break;
					}	
				}
			} 
		#endif			
	}
	fprintf(stderr,"Fragment Extraction finished.\n");

	all=malloc(sizeof(FRAGMENTS));
	for(i=0;i<MaxRes;i++)
	{
		all->res[i]= (RESIDUE *) malloc(sizeof(RESIDUE));
		for(j=0;j<MaxAtom;j++)
			all->res[i]->atom[j]=(ATOM *) malloc(sizeof(ATOM));
	}



	/*** FRAGMENT VALIDATION ****/
	if ( (Atmrec.frag[0]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}
	if ( (Atmrec.frag[1]= (FRAGMENTS *) calloc(1,sizeof(FRAGMENTS))) == NULL ) { printf("Error in malloc\n");}

	for(i=0;i<m-5;i++)
	{
		#if EXHAUSTIVE
		for(j=Worst[i];j<100;j++)
		{	
			for(new_frag= Frag_lib[i][j];new_frag!=NULL;new_frag=new_frag->next)
			{
				top=0;
				sscanf(new_frag->line,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d",Header,&Chain,&start1,&k,DB_Seq,&c,&ramach_score,&seq_score,&length,&start2,&resolution,&ss_score);
				//printf("%d\n",length);
				//printf("%s\t%c\t%d\t%d\t%s\t%c\n",Header,Chain,start1,length,DB_Seq,c);
				Atmrec.frag[0]->start_res	= start1;
				Atmrec.frag[1]->start_res	= start2;
				sprintf(Atmrec.frag[0]->fname,"%s/%c%c%c%c.pdb",getenv("PDB"),tolower(Header[0]),tolower(Header[1]),tolower(Header[2]),tolower(Header[3]));
				top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,Chain,DB_Seq);
				
				if(VALIDATE)
				{
					sprintf(Atmrec.frag[1]->fname,"./%s.pdb",argv[1]);
					top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,'A',NULL);
				}
				if(top<1+VALIDATE) /* If, by any reason, was not capable of finding the second fragment! */
					fprintf(stderr,"Could not open the fragment starting at position %d for %s_%c.\n",start1,Header,Chain);
				else
				{
					if(VALIDATE)
					{
						if( Atmrec.frag[0]->natom == Atmrec.frag[1]->natom  &&  Atmrec.frag[1]->natom > 0  &&  Atmrec.frag[0]->natom > 0 && Atmrec.frag[0]->numres > 0 && Atmrec.frag[1]->numres > 0 && Atmrec.frag[0]->numres == Atmrec.frag[1]->numres )
						{
							rmsd =  super4( Atmrec.frag[0], Atmrec.frag[1],length);
							if(rmsd > 10.0 || rmsd < 0.0) continue;
							sprintf(AUX,"%.2lf",rmsd);
							new_frag->line[strlen(new_frag->line)-1]='\t';
							strcat(new_frag->line,AUX);	
						}
					}
					if( coordcol2(all,Header,Chain,DB_Seq,&start_res) ) 
						continue;
					new_score=0.0;
					for(k=0, k2=start_res ; k2 < start_res+length ; k++, k2++)
					{
						if(k2<1 || k2 >= all->numres-1)
							continue;
						Phi[k] = Calc_dih( all->res[k2-1] , all->res[k2] , all->res[k2+1], &Psi[k]);
						new_score+= sqrt(pow(Phi[k]-Angles[start2+k][0],2)/Angles[start2+k][2] + pow(Psi[k]-Angles[start2+k][1],2)/Angles[start2+k][3]);
					}
					if(!VALIDATE)
					{
							new_frag->line[strlen(new_frag->line)-1]='\t';
							printf("%s%.2lf\n",new_frag->line,new_score);	
					}
					else
						printf("%s\t%.2lf\n",new_frag->line,new_score);	
				}
			}
		}
		#endif

		#if RANDOM
		for(j=99;j>=0;j--)
		{
			for(new_frag= Frag_lib_rnd[i][j];new_frag!=NULL;new_frag=new_frag->next)
			{
				top=0;
				sscanf(new_frag->line,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d",Header,&Chain,&start1,&k,DB_Seq,&c,&ramach_score,&seq_score,&length,&start2,&resolution,&ss_score);
				Atmrec.frag[0]->start_res	= start1;
				Atmrec.frag[1]->start_res	= start2;

				sprintf(Atmrec.frag[0]->fname,"%s/%c%c%c%c.pdb",getenv("PDB"),tolower(Header[0]),tolower(Header[1]),tolower(Header[2]),tolower(Header[3]));
				top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,Chain,DB_Seq);
				if(VALIDATE)
				{				
	                                sprintf(Atmrec.frag[1]->fname,"./%s.pdb",argv[1]);
					top = coordcol(Atmrec.frag[top]->fname,top,length,Atmrec.frag[top]->start_res,'A',NULL);
				}

				if(top<1+VALIDATE) /* If, by any reason, was not capable of finding the second fragment! */
					fprintf(stderr,"Could not open the fragment starting at position %d for %s_%c.\n",start1,Header,Chain);
				else
				{
					if(VALIDATE)
					{
						if( Atmrec.frag[0]->natom == Atmrec.frag[1]->natom  &&  Atmrec.frag[1]->natom > 0  &&  Atmrec.frag[0]->natom > 0 && Atmrec.frag[0]->numres > 0 && Atmrec.frag[1]->numres > 0 && Atmrec.frag[0]->numres == Atmrec.frag[1]->numres )
						{
							rmsd =  super4( Atmrec.frag[0], Atmrec.frag[1],length);
							if(rmsd > 10.0 || rmsd < 0.0) continue;
							sprintf(AUX,"%.2lf",rmsd);
							new_frag->line[strlen(new_frag->line)-1]='\t';
							strcat(new_frag->line,AUX);	
						}					
					}	
					
					if( coordcol2(all,Header,Chain,DB_Seq,&start_res) ) 
						continue;
					new_score=0.0;
					for(k=0, k2=start_res ; k2 < start_res+length && k2 < all->numres; k++, k2++)
					{
						if(k2<1)
							continue;
						Phi[k] = Calc_dih( all->res[k2-1] , all->res[k2] , all->res[k2+1], &Psi[k]);
						new_score+= sqrt(pow(Phi[k]-Angles[start2+k][0],2)/Angles[start2+k][2] + pow(Psi[k]-Angles[start2+k][1],2)/Angles[start2+k][3]);
					}
					
					if(!VALIDATE)
					{
							new_frag->line[strlen(new_frag->line)-1]='\t';
							printf("%s%.2lf\n",new_frag->line,new_score);	
					}
					else
						printf("%s\t%.2lf\n",new_frag->line,new_score);	
				}
			}
		}
		#endif
	}

	free(Atmrec.frag[0]);
	free(Atmrec.frag[1]);

	free(Total);
	free(Total_rnd);
	free(Worst);
	free(Worst_rnd);
	fclose(input_fasta);
	fclose(input_ss);
	fclose(input_pdb);
	fclose(input_phipsi);
	fclose(blossum_file);
	return 0;
}


