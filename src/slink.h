#define FILEL 50
#define MaxFrag 10000


#define TRUE 1
#define FALSE 0


typedef struct {
char Ffile[FILEL];
char Pfile[FILEL];
char Outfil[FILEL];
} Car;


typedef struct {
char Line[MaxFrag][100];
char fname[MaxFrag][FILEL];
int start_res[MaxFrag];
int Topf;
int length;
} Data;


#define  MaxRes  1000 /* Maximum  no. of residues total */
#define  MaxAtom 6    /* Max number of atoms per residue */
#define  MAX_FIT_COORDS 200


typedef struct
{
 char atomname[4];
 float x;
 float y;
 float z;
} ATOM;

/* Each RESIDUE  has numatom atoms,  residue name and number */
typedef struct
{
  ATOM  * atom[MaxAtom];
  char Resname[4];
  short int Resno;
  short int Resindex;
  short int numatom;
  char type;    /* type of actual/predicted secondary structure */
} RESIDUE;


/* Each FRAGMENT has numres RESIDUES */
typedef struct
{
  RESIDUE * res[MaxRes];
  int numres;
  int natom;
  char fname[FILEL];
  char Header[5];
  char Chain;
  char ALN_Seq[25];
  char aux;
  int match_score;
  int seq_score; 
  int length; 
  int start2; 
  double resolution;
  double rmsd;
  int ss_score;
  int start_res; 
  int end_res;
  int printed;
} FRAGMENTS;

typedef struct
{
  FRAGMENTS * frag[MaxFrag];
}BITS ;

