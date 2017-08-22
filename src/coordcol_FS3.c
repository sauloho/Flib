#include "slink.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*C Deane 27/10/04 Collect Coordinates*/

char Dic1[24][4]={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
				  "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","ASX","GLX",
				  "BUM","BUL"	};      /* A 3-letter code residue dictionary.              */


char Dic2[24]={'A','R','N','D','C','Q','E','G','H','I','L',
               'K','M','F','P','S','T','W','Y','V','B','Z',
               'X','*'};      /* A one-letter code residue dictionary.              */

extern BITS Atmrec;

char three_to_one(char three[4])
{
	int i;
	for(i=0;i<22;i++)
	{
		if(!strcmp(three,Dic1[i]))
			return Dic2[i];
	}
	return '*';
}

void extract_frags(char bum[FILEL], int length, int start_res,char ChainQ,char* Seq)
{
	int a;
	char line[150],chain;
	char* string_ptr;
	float xcoord,ycoord,zcoord;
	char name[9];
	char type[5];
	char restype[4];
	char Sequence[2000];
	FILE *IN,*temp;
	int count,count_aux,old_count_aux;
	int position;
	char * pch;

	count=-1;
	old_count_aux = -1;

	/*******************************************************/

	temp=fopen("temp.pdb","w");

	/* First look for the position on the sequence of the file */
	if(Seq!=NULL)
	{	
		IN=fopen(bum,"r");
		if (IN == NULL)
		{
	    	    (void)printf("unable to open file %s in coordcol\n",bum);  /* error message if file unavailable */
 	   	     return ;
    		    /*exit(8); Exit to operating system */
    		}

		while(!feof(IN))
		{
			if(count>=2000)
				break;
			string_ptr = fgets(line,sizeof(line),IN);
			line[strlen(line)-1] = '\0';
			sscanf(line,"%4s%*9c%3s%*1c%3s%*1c%1c%*1c%3d%*5c%9f%9f%9f\n",type,name,restype,&chain,&count_aux,&xcoord,&ycoord,&zcoord);
			if( ( !strcmp(type,"ATOM") ) && chain == ChainQ )
			{
				if( ( !strcmp(name,"CA") ) && old_count_aux != count_aux)
				{
					count++;
					Sequence[count]=three_to_one(restype);
				}
			}
		}
		fclose(IN);
	
		pch = strstr (Sequence,Seq);
		if(pch!=NULL)
		{
			position = start_res;
			start_res = pch - Sequence;
		}
		else
			return ;
	}
	
	/*******************************************************/

	count=-1;
	old_count_aux = -1;

	IN=fopen(bum,"r");
	if (IN == NULL)
	{
        	(void)printf("unable to open file %s in coordcol\n",bum);  /* error message if file unavailable */
	         return ;
        	/*exit(8); Exit to operating system */
	}

	while(!feof(IN))
	{
		string_ptr = fgets(line,sizeof(line),IN);
		line[strlen(line)-1] = '\0';
		//fprintf(stderr,"%s\n",line);
		sscanf(line,"%4s%*9c%3s%*1c%3s%*1c%1c%*1c%3d%*5c%9f%9f%9f\n",type,name,restype,&chain,&count_aux,&xcoord,&ycoord,&zcoord);
		//fprintf(stderr,"%4s %3s %3s %1c %3d %9f %9f %9f\n",type,name,restype,chain,count_aux,xcoord,ycoord,zcoord);		
		if( ( !strcmp(type,"ATOM") ) && chain == ChainQ )
		{
		    if( old_count_aux != count_aux)
		    {	
	  	    	    old_count_aux = count_aux;	/* This is a new residue! YAY! */
			    count++;                 	/* Increase the CA counter for the chain... */		
		    }
		    if(count_aux >= start_res && count_aux < start_res+length)
				fprintf(temp,"%s\n",line);	
		}
 	}/*while*/
	fclose(temp);
	fclose(IN);
}
