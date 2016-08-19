#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>     

int main(int argc, char *argv[])
{
	int length,start1,match_score,seq_score,start2,ss_score,total,old=-1;
	char Header[82],Chain,type,ALN_Seq[30],c;
	char Fasta_Seq[1000];
	char AUX[82];
	float rmsd,score;
	double resolution;
	FILE *library;
	FILE *input_fasta;
	FILE *out_500,*out_20;

	if(argc <= 2 )
	{
		printf("USAGE: %s PDB_ID FLIB_FILE.lib3000\n",argv[0]);
		return 0;
	}

  	/*** FILE HANDLING: ***/
	library = fopen(argv[2],"r");
	if(library==NULL) { fprintf(stderr,"Library file: %s not found!\n",argv[2]); 	return 0; 	}

	strcpy(AUX,argv[1]);
	input_fasta = fopen(strcat(AUX,".fasta.txt"),"r");
	if(input_fasta==NULL) { fprintf(stderr,"Fasta file: %s not found!\n",AUX); 	return 0;	}

	/* READ INPUT FASTA SEQUENCE */
    /* Remove Header from FASTA file */
    for(c = fgetc(input_fasta); c!='\n' ; c=fgetc(input_fasta));
    fscanf(input_fasta,"%s",Fasta_Seq);

	while(fscanf(input_fasta,"%s",AUX)!=EOF) 	
		strcat(Fasta_Seq,AUX);


	strcpy(AUX,argv[1]);
	out_20 = fopen(strcat(AUX,".lib20"),"w");

    strcpy(AUX,argv[1]);
    out_500 = fopen(strcat(AUX,".lib500"),"w");

	for(total=0;fscanf(library,"%s\t%c\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t%lf\t%d\t%f\n",Header,&Chain,&start1,&length,ALN_Seq,&type,&match_score,&seq_score,&length,&start2,&resolution,&ss_score,&rmsd)!=EOF;total++)
	{
		if(old != start2) /* Fragments for a new position... */
		{
			total=0;
			old = start2;
		}
        if(total<500)
		        fprintf(out_500,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,rmsd);
		if(total<20)
			fprintf(out_20,"%s\t%c\t%3d\t%3d\t%s\t%c\t%3d\t%3d\t%3d\t%3d\t%.2lf\t%d\t%.2f\n",Header,Chain,start1,length,ALN_Seq,type,match_score,seq_score,length,start2,resolution,ss_score,rmsd);
	}
	fclose(library);
	fclose(input_fasta);
	fclose(out_20);
	fclose(out_500);
    
	return 0;
}
