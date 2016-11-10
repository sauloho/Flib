/*
Originally written in Fortran by B.S.Neela and H.A.Nagarajaram

For the method the reference you should quote is:


SIMON K. KEARSELY ACTA.CRYST.(1989), A45, 208-210.


Obviously the two molecules A and B should have the same
number of atoms


Conversion to C restart arrays at zero by Charlotte Deane */


/*editing for use in FRAG SUB C Deane
28/10/04
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jacobi.h"
#include "slink.h"


extern BITS Atmrec;

int nrot;

void eig(double a[4][4], double d[4], double v[4][4]);


float super4(FRAGMENTS *One, FRAGMENTS *Two, int nres)
{
float ormsd;
int n;
double A[4][4];
double D[4];
double V[4][4];
double P[4][4];
double XORT[MAX_FIT_COORDS],YORT[MAX_FIT_COORDS],ZORT[MAX_FIT_COORDS];
double XXORT[MAX_FIT_COORDS],YYORT[MAX_FIT_COORDS],ZZORT[MAX_FIT_COORDS];
double XR[MAX_FIT_COORDS],YR[MAX_FIT_COORDS],ZR[MAX_FIT_COORDS];
double R[3][3];
double X1_OUT[3][MAX_FIT_COORDS],X2_OUT[3][MAX_FIT_COORDS];
int i,j,k,a,b;
int natom;
double store[500][3];



double SUMXORT=0.0;
double SUMYORT=0.0;
double SUMZORT=0.0;
double SUMXXORT=0.0;
double SUMYYORT=0.0;
double SUMZZORT=0.0;

double CXX,CYY ,CZZ,CX,CY,CZ;
double XM,YM,ZM,XP,YP,ZP;
double Q1,Q2,Q3,Q4;
k=0;
for(i=0;i<nres;i++){
  for(j=0;j<One->res[i]->numatom;j++){
    if(!strcmp(One->res[i]->atom[j]->atomname,"CA")){
      XORT[k] = One->res[i]->atom[j]->x;
      YORT[k] = One->res[i]->atom[j]->y;
      ZORT[k] = One->res[i]->atom[j]->z;
      /*printf("%f %f %f x,y,z %d k\n",XORT[k], YORT[k],ZORT[k],k);*/
      k++;
    }
  }
}
natom =k;

k=0;
for(i=0;i<nres;i++){
  for(j=0;j<Two->res[i]->numatom;j++){
    if(!strcmp(Two->res[i]->atom[j]->atomname,"CA")){
      XXORT[k] = Two->res[i]->atom[j]->x;
      YYORT[k] = Two->res[i]->atom[j]->y;
      ZZORT[k] = Two->res[i]->atom[j]->z;
      /*printf("%f %f %f x,y,z %d k\n",XXORT[k], YYORT[k],ZZORT[k],k);*/
      k++;
    }
  }
}

if(k!=natom){fprintf(stderr,"different number of atoms in super4???\n"); return -1.0;}
/*
for(i=0;i<natom;i++){
XORT[i] = X1_IN[0][i];
YORT[i] = X1_IN[1][i];
ZORT[i] = X1_IN[2][i];

XXORT[i] = X2_IN[0][i];
YYORT[i] = X2_IN[1][i];
ZZORT[i] = X2_IN[2][i];
}
*/

for(i=0;i<natom;i++){
SUMXORT=SUMXORT+XORT[i];
SUMYORT=SUMYORT+YORT[i];
SUMZORT=SUMZORT+ZORT[i];
SUMXXORT=SUMXXORT+XXORT[i];
SUMYYORT=SUMYYORT+YYORT[i];
SUMZZORT=SUMZZORT+ZZORT[i];
}


CXX = SUMXXORT/(natom*1.0);
CYY = SUMYYORT/(natom*1.0);
CZZ = SUMZZORT/(natom*1.0);
CX = SUMXORT/(natom*1.0);
CY = SUMYORT/(natom*1.0);
CZ = SUMZORT/(natom*1.0);



for(i=0;i<4;i++){
  for(j=0;j<4;j++){
	P[i][j]=0.0;
  }
}
for (i=0;i<natom;i++){
  XORT[i]=XORT[i]-CX;
  YORT[i]=YORT[i]-CY;
  ZORT[i]=ZORT[i]-CZ;
  XXORT[i]=XXORT[i]-CXX;
  YYORT[i]=YYORT[i]-CYY;
  ZZORT[i]=ZZORT[i]-CZZ;

  XM=XORT[i]-XXORT[i];
  YM=YORT[i]-YYORT[i];
  ZM=ZORT[i]-ZZORT[i];
  XP=XORT[i]+XXORT[i];
  YP=YORT[i]+YYORT[i];
  ZP=ZORT[i]+ZZORT[i];



/*
     LEAST SQUARES FIT ALGEBRA GIVEN IN ACTA.CRYST.(1989)
     A45 208-210. BY SIMON K. KEARSELY

*/

 P[0][0]=P[0][0]+(XM*XM+YM*YM+ZM*ZM);
 P[0][1]=P[0][1]+(YP*ZM-YM*ZP);
 P[0][2]=P[0][2]+(XM*ZP-XP*ZM);
 P[0][3]=P[0][3]+(XP*YM-XM*YP);
 P[1][0]=P[0][1];
 P[1][1]=P[1][1]+(YP*YP+ZP*ZP+XM*XM);
 P[1][2]=P[1][2]+(XM*YM-XP*YP);
 P[1][3]=P[1][3]+(XM*ZM-XP*ZP);
 P[2][0]=P[0][2];
 P[2][1]=P[1][2];
 P[2][2]=P[2][2]+(XP*XP+ZP*ZP+YM*YM);
 P[2][3]=P[2][3]+(YM*ZM-YP*ZP);
 P[3][0]=P[0][3];
 P[3][1]=P[1][3];
 P[3][2]=P[2][3];
 P[3][3]=P[3][3]+(XP*XP+YP*YP+ZM*ZM);
} 



eig(P,D,V);


ormsd = sqrt(fabs(D[3]/natom));


/*printf("%f rmsd\n",ormsd);*/



/*if(ormsd < 1){*/

 Q1=V[0][3];
 Q2=V[1][3];
 Q3=V[2][3];
 Q4=V[3][3];

/*
     ROTATION MATRIX ALGEBRA GIVEN IN ACTA.CRYST.(1989) A45 208-210
*/
 


R[0][0]=Q1*Q1+Q2*Q2-Q3*Q3-Q4*Q4;
R[0][1]=2*(Q2*Q3+Q1*Q4);
R[0][2]=2*(Q2*Q4-Q1*Q3);
R[1][0]=2*(Q2*Q3-Q1*Q4);
R[1][1]=Q1*Q1+Q3*Q3-Q2*Q2-Q4*Q4;
R[1][2]=2*(Q3*Q4+Q1*Q2);
R[2][0]=2*(Q2*Q4+Q1*Q3);
R[2][1]=2*(Q3*Q4-Q1*Q2);
R[2][2]=Q1*Q1+Q4*Q4-Q2*Q2-Q3*Q3;


/*put in calculating the fitted fragment just so I can look should I care at 
a later stage*/
/*
for(i=0;i<nres;i++){
  for(j=0;j<Two->res[i]->numatom;j++){
	store[i][0] = Two->res[i]->atom[j]->x - CXX;
	store[i][1] = Two->res[i]->atom[j]->y - CYY;
	store[i][2] = Two->res[i]->atom[j]->z - CZZ;

	  Two->res[i]->atom[j]->x=R[0][0]*store[i][0]+R[0][1]*store[i][1]+R[0][2]*store[i][2];
  	  Two->res[i]->atom[j]->x=   Two->res[i]->atom[j]->x +CX;
  	  Two->res[i]->atom[j]->y=R[1][0]*store[i][0]+R[1][1]*store[i][1]+R[1][2]*store[i][2];
  	  Two->res[i]->atom[j]->y=   Two->res[i]->atom[j]->y +CY;
  	  Two->res[i]->atom[j]->z=R[2][0]*store[i][0]+R[2][1]*store[i][1]+R[2][2]*store[i][2];
  	  Two->res[i]->atom[j]->z=   Two->res[i]->atom[j]->z +CZ;

  }
}

*/


/*}  if ORMSD*/

/*printf("%f rmsd\n",ormsd);*/
return(ormsd);
}


/*--------------------------------------------------------------
    
     THE FOLLOWING THREE SUBROUTINES ARE FROM THE BOOK
                  NUMERICAL RECIPES IN C ON LINE
    pg463 -469 

  
-----------------------------------------------------------*/



void eig(double A[4][4], double D[4], double V[4][4])
{
jacobi(A,4,D,V,nrot);
}
 


						
