#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define MAXLEN 10000

int main(int argc, char *argv[])
{
	int i,j,k,cov,cov2,h,b,l,good_h,good_b,good_l,end;
	int total;
	int Covering[MAXLEN],Covering2[MAXLEN],numgood=0;
	int Numgood[4],Total[4];
	char AUX[82];
	char Fasta_Seq[MAXLEN];		  /* The protein's fasta sequence.					    */
	int length,start1,match_score,seq_score,start2,ss_score;
	char Header[82],Chain,type,ALN_Seq[200];
	float rmsd,avg_length=0.0,torsion;
	double resolution, CUTOFF;
  	char c,DSSP_Seq[1000],SS[1000];		  /* The query protein's secondary structure			*/	
	FILE *input_fasta;
	FILE *library;
	FILE *input_dssp;

	if(argc<3)
	{
		printf("USAGE: %s FRAG_LIB PDB_ID CUTOFF\n",argv[0]);
		return 0;
	}
  	/*** FILE HANDLING: ***/
	library = fopen(argv[1],"r");
	CUTOFF = atof(argv[3]);	
	for(i=0;i<4;i++)
	{
		Numgood[i]=0;
		Total[i]=0;
	}

	/* READ INPUT SECONDARY STRUCTURE SEQUENCE */		
	strcpy(AUX,argv[2]);
	input_dssp = fopen(strcat(AUX,".dssp_psi"),"r");
	fscanf(input_dssp,"%s",DSSP_Seq);

	//printf("%s   %d\n",DSSP_Seq,strlen(DSSP_Seq));
	
	l=0; h=0; b=0; cov2=0;
	for(i=0;i<strlen(DSSP_Seq);i++)
	{
		switch(DSSP_Seq[i])
		{
			case 'H':
				h++;
				SS[i]='H';			
				break; 
			case 'E':
				b++;
				SS[i]='B';
				break;
			case 'B':
				b++;
				SS[i]='B';
				break;
			default:
				l++;
				SS[i]='L';
				break;
		}
	}	
	//printf("%s\n",SS);

	/* READ INPUT FASTA SEQUENCE */
	input_fasta = fopen(strcat(argv[2],".fasta.txt"),"r");
	fscanf(input_fasta,"%s",Header);
	for(c=fgetc(input_fasta);c!='\n';c=fgetc(input_fasta));
	fscanf(input_fasta,"%s",Fasta_Seq);
	while(fscanf(input_fasta,"%s",AUX)!=EOF) 	
		strcat(Fasta_Seq,AUX);
	
	total = strlen(Fasta_Seq);
	for(i=0;i<total;i++)
	{
		Covering[i]=0;
		Covering2[i]=0;
	}
	numgood=0;


	for(j=0;fscanf(library,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\t%f\n",Header,&Chain,&start1,&end,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&torsion,&rmsd)!=EOF;j++)
	{
		if(rmsd < 0.0)
		{
			continue;
		}
		avg_length+=length;
		for(k=start2; k < start2+length && k<total;k++)
			Covering2[k]=1;
		/* Get actual coverage and precision */
		if(rmsd < CUTOFF)
		{
			numgood++;
			for(k=start2 ; k < start2 + length && k < total; k++)
				Covering[k]=1;
			switch (type)
            {
                 case 'H':
			        Numgood[0]++;
                    break;
                 case 'B':
                    Numgood[1]++;
					break;
                 case 'E':
                    Numgood[1]++;
					break;
                 case 'O':
                    Numgood[2]++;
					break;
                 case 'L':
                    Numgood[3]++;
					break;
                 case ' ':
                    Numgood[3]++;
					break;
	             default:	
                    Numgood[3]++;
					break;
            }
		}
		
		switch (type)
        {
             case 'H':
		        Total[0]++;
                break;
             case 'B':
                Total[1]++;
				break;
             case 'E':
                Total[1]++;
				break;
             case 'O':
                Total[2]++;
				break;
             case 'L':
                Total[3]++;
				break;
			case ' ':
                Total[3]++;
				break;
             default:	
				Total[3]++;
				break;
         }
	}

	
	good_l=0; good_h=0; good_b=0; 
	for(k=0,cov=0,cov2=0;k<total;k++)
	{
		if(Covering2[k]) cov2++;
		if(Covering[k])
		{
			switch(SS[k])
			{
		             case 'H':
				good_h++;
	                	break;
		             case 'B':
				good_b++;
	                	break;
	        	     case 'L':
				good_l++;
	        	        break;
			}
			cov++;
		}
	}
	start2++;

	printf("%.2f\t%.2f\t%.2f",(float)j/start2,(float)numgood/j,(float)cov/cov2);
	for(i=0;i<4;i++)
		printf("\t%.3f",(float)Numgood[i]/(Total[i]));
	printf("\t%.3f\t%.3f\t%.3f\t%.3f\n",(float)good_h/h,(float)good_b/b,(float)good_l/l,avg_length/(float)j);
	return 0;
}
