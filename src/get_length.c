#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc, char *argv[])
{
	FILE *input_fasta;
	char AUX[202];
	char Header[202];              /* Information extracted from the header of the fasta */
	char Fasta_Seq[5000];		  /* The protein's fasta sequence.					    */

	if(!argv[1])
	{
		printf("USAGE: %s target_fasta_file.fasta\n",argv[0]);
		return 0;
	}
	input_fasta = fopen(argv[1],"r");
	/* READ INPUT FASTA SEQUENCE */
	fscanf(input_fasta,"%s",Header);
	fscanf(input_fasta,"%s",Fasta_Seq);
	while(fscanf(input_fasta,"%s",AUX)!=EOF) 	
		strcat(Fasta_Seq,AUX);
	printf("%d\n",strlen(Fasta_Seq));

	return 0;
}
