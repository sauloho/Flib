#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "slink.h"

#define MAXLEN 10000
#define MaxAtom 6    /* Max number of atoms per residue */

static char Dic1[24][4]={"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU",
				  "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","ASX","GLX",
				  "BUM","BUL"	};      /* A 3-letter code residue dictionary.              */


static char Dic2[24]={'A','R','N','D','C','Q','E','G','H','I','L',
               'K','M','F','P','S','T','W','Y','V','B','Z',
               'X','*'};      /* A one-letter code residue dictionary.              */

static char three_to_one(char three[4])
{
	int i;
	for(i=0;i<22;i++)
	{
		if(!strcmp(three,Dic1[i]))
			return Dic2[i];
	}
	return '*';
}

short int ResIndex(char  * Resname){
             if( strncmp(Resname,"GLY",3) == 0 ) return 0;
             if( strncmp(Resname,"ALA",3) == 0 ) return 1;
             if( strncmp(Resname,"VAL",3) == 0 ) return 2;
             if( strncmp(Resname,"LEU",3) == 0 ) return 3;
             if( strncmp(Resname,"ILE",3) == 0 ) return 4;
             if( strncmp(Resname,"PRO",3) == 0 ) return 5;
             if( strncmp(Resname,"ASP",3) == 0 ) return 6;
             if( strncmp(Resname,"GLU",3) == 0 ) return 7;
             if( strncmp(Resname,"ASN",3) == 0 ) return 8;
             if( strncmp(Resname,"GLN",3) == 0 ) return 9;
             if( strncmp(Resname,"LYS",3) == 0 ) return 10;
             if( strncmp(Resname,"ARG",3) == 0 ) return 11;
             if( strncmp(Resname,"SER",3) == 0 ) return 12;
             if( strncmp(Resname,"THR",3) == 0 ) return 13;
             if( strncmp(Resname,"MET",3) == 0 ) return 14;
             if( strncmp(Resname,"CYS",3) == 0 ) return 15;
             if( strncmp(Resname,"TYR",3) == 0 ) return 16;
             if( strncmp(Resname,"TRP",3) == 0 ) return 17;
             if( strncmp(Resname,"HIS",3) == 0 ) return 18;
             if( strncmp(Resname,"PHE",3) == 0 ) return 19;
             if( strncmp(Resname,"UNK",3) == 0 ) return -1;
             if( strncmp(Resname,"PCA",3) == 0 ) return -2;
             if( strncmp(Resname,"ACE",3) == 0 ) return -3;
             if( strncmp(Resname,"HOH",3) == 0 ) return -4;
/*             printf("ERROR IN ResIndex - unknown 
Resnameidue=%s\n",Resname);*/
             return -5;
}


double dihed(double list[3][4])
{
	int j;
	double co1[3];
	double co2[3];
	double co3[3];
	double co4[3];
	/*need to use a set of vectors to calculate the dihedral*/
	double v_12[3]; /*vector from 1 to 2 */
	double v_23[3]; /*vector from 2 to 3 */
	double v_34[3]; /*vector from 3 to 4 */
	double n_1[3]; 	/*normal vector of the plane 1-2-3*/
	double n_2[3]; 	/*normal vector of the plane 2-3-4*/
	double abs_n1;  /*length of the vectors*/
	double abs_n2;
	double abs_n1_n2;
	double n1_n2; 	/*scalar product (dot product) of n_1 and n_2*/
	double ss; 		/*scalar product of v_12 and n_2 to determine the sign of phi */
	double dihedral; /*the dihedral angle*/

	for(j=0;j<3;j++)
	{
		co1[j]=list[j][0];
		co2[j]=list[j][1];
		co3[j]=list[j][2];
		co4[j]=list[j][3];
	}

	for(j=0;j<3;j++)
	{
		v_12[j] = co2[j] - co1[j];
		v_23[j] = co3[j] - co2[j];
		v_34[j] = co4[j] - co3[j];
	}

	n_1[0] = v_12[1] * v_23[2] - v_12[2] * v_23[1];
	n_1[1] = v_12[2] * v_23[0] - v_12[0] * v_23[2];
	n_1[2] = v_12[0] * v_23[1] - v_12[1] * v_23[0];

	abs_n1 = n_1[0] * n_1[0] + n_1[1] * n_1[1] + n_1[2] *  n_1[2];
	n_2[0] = v_23[1] * v_34[2] - v_23[2] * v_34[1];
	n_2[1] = v_23[2] * v_34[0] - v_23[0] * v_34[2];
	n_2[2] = v_23[0] * v_34[1] - v_23[1] * v_34[0];

	abs_n2 = n_2[0] * n_2[0] + n_2[1] * n_2[1] + n_2[2] *  n_2[2];
	abs_n1_n2 = abs_n1 * abs_n2;

	if(abs_n1_n2 < 0 || abs_n1_n2 != abs_n1_n2 )
	{
		fprintf(stderr,"warning: dihedral angle not defined \n");
		dihedral =0;
		return 0;
	}

	n1_n2 = n_1[0] * n_2[0] + n_1[1] * n_2[1] + n_1[2] * n_2[2];
	ss = v_12[0] * n_2[0] + v_12[1] * n_2[1] +v_12[2] * n_2[2];

	if(ss >= 0.)
  		dihedral = acos(n1_n2/sqrt(abs_n1_n2));
	else
		dihedral = -acos(n1_n2/sqrt(abs_n1_n2));

	dihedral = (180.0/3.14159) * dihedral;

	if (dihedral!=dihedral)
		return 0;
	return dihedral;
}



float Calc_dih(RESIDUE * OneP,RESIDUE * One,RESIDUE * OneA, float * temp)
{
	double List_phi[3][4];
	double List_psi[3][4];
	int a;
	float phi;

	/*printf("%d OneP->numatom\n",OneP->numatom);*/
	for(a=0;a<OneP->numatom;a++)
	{
		if(!strcmp(OneP->atom[a]->atomname,"C"))
		{
			List_phi[0][0] = OneP->atom[a]->x;
			List_phi[1][0] = OneP->atom[a]->y;
			List_phi[2][0] = OneP->atom[a]->z;
			break;
		}
	}

	for(a=0;a<One->numatom;a++)
	{
		if(!strcmp(One->atom[a]->atomname,"N"))
		{
			List_phi[0][1] = One->atom[a]->x;
			List_phi[1][1] = One->atom[a]->y;
    			List_phi[2][1] = One->atom[a]->z;
	    		List_psi[0][0] = One->atom[a]->x;
    			List_psi[1][0] = One->atom[a]->y;
    			List_psi[2][0] = One->atom[a]->z;
  		}
		else if(!strcmp(One->atom[a]->atomname,"CA"))
		{
	     		List_phi[0][2] = One->atom[a]->x;
			List_phi[1][2] = One->atom[a]->y;
			List_phi[2][2] = One->atom[a]->z;
			List_psi[0][1] = One->atom[a]->x;
			List_psi[1][1] = One->atom[a]->y;
			List_psi[2][1] = One->atom[a]->z;
		}
		else if(!strcmp(One->atom[a]->atomname,"C"))
		{
			List_phi[0][3] = One->atom[a]->x;
			List_phi[1][3] = One->atom[a]->y;
			List_phi[2][3] = One->atom[a]->z;
			List_psi[0][2] = One->atom[a]->x;
			List_psi[1][2] = One->atom[a]->y;
			List_psi[2][2] = One->atom[a]->z;
		}
	}

	for(a=0;a<OneA->numatom;a++)
	{
		if(!strcmp(OneA->atom[a]->atomname,"N"))
		{
		    List_psi[0][3] = OneA->atom[a]->x;
		    List_psi[1][3] = OneA->atom[a]->y;
		    List_psi[2][3] = OneA->atom[a]->z;
			break;
		 }	
	}

	*temp = dihed(List_psi);
	phi = dihed(List_phi);

	return phi;
}

int coordcol2(FRAGMENTS * Piece,char PDB[82],char Chain, char* Seq,int *start)
{
	int b;
	int i;
	char line[500];
	char Sequence[MaxRes];
	int start_res,old_count_aux;
	char * pch;	
	char* string_ptr;
	float xcoord,ycoord,zcoord;
	char name[4];
	char type[5];
	char restype[4];
	char chain;
	FILE *IN;
	int count,count_aux;
	int oldcount;
	int natom;
	int w=-1;
	char filename[120];

	sprintf(filename,"%s/",getenv("PDB"));
	for (i=0 ;i<4; i++) PDB[i] = tolower(PDB[i]);
	strcat(filename,PDB);
	strcat(filename,".pdb");
	IN = fopen(filename,"r");
	if(IN==NULL)
	{
		fprintf(stderr,"Could not find the file %s\n",filename);
		return 1;	
	}

	count=-1;
	old_count_aux = -999999;
	oldcount=-999999;
	b=0;
	natom=0;
	w=-1;

	/* First look for the position on the sequence of the file */
	if(Seq!=NULL)
	{	
		while(!feof(IN))
		{
			string_ptr = fgets(line,sizeof(line),IN);
			sscanf(line,"%4s%*9c%3s%*1c%3s%*1c%1c%*1c%3d%*5c%9f%9f%9f\n",type,name,restype,&chain,&count_aux,&xcoord,&ycoord,&zcoord);
			if( ( !strcmp(type,"ATOM") ) && chain == Chain )
			{
				if((!strcmp(name,"CA")) || (!strcmp(name,"C")) || (!strcmp(name,"N")) || (!strcmp(name,"O")) || (!strcmp(name,"CB")))
				{
					if( old_count_aux != count_aux )
					{
						if(count >= 0 && count < MaxRes)
						    b >=5 ? Piece->res[count]->numatom=5 : b;

						old_count_aux = count_aux;
						count++;
						if(count >= MaxRes)
						{
							fprintf(stderr,"Warning: PDB file %s has more than %d residues! Sequence discarded.\n",PDB,MaxRes);
							fclose(IN);
							return 1;
						}
						Sequence[count]=three_to_one(restype);

          				b=0;	
	        			strcpy((Piece->res[count])->Resname,restype);
    	    			(Piece->res[count])->Resno = count_aux;
						(Piece->res[count])->Resindex = ResIndex((Piece->res[count])->Resname);
		        	}

			    	(Piece->res[count])->atom[b]->x = xcoord;	(Piece->res[count])->atom[b]->y = ycoord;	(Piece->res[count])->atom[b]->z = zcoord;
				    strcpy((Piece->res[count])->atom[b]->atomname,name);

	        		if(b < 5)
					{      
	        			b++;
		        		natom++;
	        		}
        		}	
			}
		}
		fclose(IN);
		Sequence[count+1]='\0';
		pch = strstr (Sequence,Seq);
		if(pch != NULL)
		{
			start_res = pch - Sequence;
			*start = start_res;
		}
		else
		{
			fprintf(stderr,"Could not find sequence %s on %s\n",Seq,PDB);
			return 1;
		}
	}

	Piece->numres = count+1;
	Piece->natom = natom;
	return 0;
}

