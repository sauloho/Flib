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

int coordcol(char bum[FILEL],int u, int length, int start_res,char ChainQ,char* Seq)
{
	int a;
	char line[150],chain;
	char* string_ptr;
	float xcoord,ycoord,zcoord;
	char name[9];
	char type[5];
	char restype[4];
	char Sequence[2000];
	FILE *IN;
	int count,count_aux,old_count_aux;
	int natom,position;
	char * pch;

	a=0;
	natom=0;
	count=-1;
	old_count_aux = -1;

	/*******************************************************/


	/* First look for the position on the sequence of the file */
	if(Seq!=NULL)
	{	
		IN=fopen(bum,"r");
		if (IN == NULL)
	    	{
    		    (void)fprintf(stderr,"unable to open file %s in coordcol\n",bum);  /* error message if file unavailable */
	    	     return u;
    		}
		while(!feof(IN))
		{
			if(count>=2000)
				break;
			string_ptr = fgets(line,sizeof(line),IN);
			line[strlen(line)-1] = '\0';
			//fprintf(stderr,"%s\n",line);
			sscanf(line,"%4s%*9c%3s%*1c%3s%*1c%1c%*1c%3d%*5c%9f%9f%9f\n",type,name,restype,&chain,&count_aux,&xcoord,&ycoord,&zcoord);
			//fprintf(stderr,"%4s %3s %3s %1c %3d %9f %9f %9f\n",type,name,restype,chain,count_aux,xcoord,ycoord,zcoord);		
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
			//fprintf(stderr,"FOUND FRAG AT POSITION %d (was supposed to be at %d)\n",start_res,position);
		}
		else
		{
			//fprintf(stderr,"DID NOT FIND FRAG!\n");
			return u;
		}
	}
	
	/*******************************************************/

	a=0;
	natom=0;
	count=-1;
	old_count_aux = -1;

	IN=fopen(bum,"r");
	if (IN == NULL)
	{
        	(void)fprintf(stderr,"unable to open file %s in coordcol\n",bum);  /* error message if file unavailable */
	         return u;
        	/*exit(8); Exit to operating system */
    	}
	//printf("%c %d %d\n",ChainQ,start_res,length);
	while(!feof(IN))
	{
		string_ptr = fgets(line,sizeof(line),IN);
		line[strlen(line)-1] = '\0';
		//fprintf(stderr,"%s\n",line);
		sscanf(line,"%4s%*9c%3s%*1c%3s%*1c%1c%*1c%3d%*5c%9f%9f%9f\n",type,name,restype,&chain,&count_aux,&xcoord,&ycoord,&zcoord);
		//fprintf(stderr,"%4s %3s %3s %1c %3d %9f %9f %9f\n",type,name,restype,chain,count_aux,xcoord,ycoord,zcoord);		
		if( ( !strcmp(type,"ATOM") ) && chain == ChainQ )
		{
			if( ( !strcmp(name,"CA") ) && old_count_aux != count_aux)
			{
				old_count_aux = count_aux;	/* This is a new residue! YAY! */
				count++;                 	/* Increase the CA counter for the chain... */	
				
			    	if(count >= start_res && count < start_res+length)
				{
					if ((Atmrec.frag[u]->res[a] = (RESIDUE *) calloc(1,sizeof(RESIDUE))) == NULL ) { fprintf(stderr,"Error in malloc\n");}
			        	if ( (Atmrec.frag[u]->res[a]->atom[0]= (ATOM *) calloc(1,sizeof(ATOM))) == NULL ) { fprintf(stderr,"Error in malloc\n");}
					(Atmrec.frag[u]->res[a])->atom[0]->x = xcoord;
				        (Atmrec.frag[u]->res[a])->atom[0]->y = ycoord;
					(Atmrec.frag[u]->res[a])->atom[0]->z = zcoord;
					strcpy((Atmrec.frag[u]->res[a])->atom[0]->atomname,name);
					strcpy((Atmrec.frag[u]->res[a])->Resname,restype);
					(Atmrec.frag[u]->res[a])->Resno = count;
					Atmrec.frag[u]->res[a]->numatom=1;
					natom++;
					a++;	
			     }
			     
 			}     	
		}
		if(a>=length)
			break;
 	}/*while*/

	fclose(IN);
	Atmrec.frag[u]->numres = a;
	//printf("%d Atmrec.frag[u] numres \n",Atmrec.frag[u]->numres);
	Atmrec.frag[u]->natom = natom;
	//printf("%d Atmrec.frag[u]->natom\n",Atmrec.frag[u]->natom);
	
	if(a!=length)
		return u;
	else 
		return(u+1);
}
