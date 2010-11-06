//------------------------------------------------------------------------------
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: show-aligns.cc
//         Date: 10 / 18 / 2002
//
//        Usage: show-aligns [options] <deltafile>
//               Try 'show-aligns -h' for more information
//
//  Description: For use in conjunction with the MUMmer package.
//              "show-aligns" displays human readable information from the
//             .delta output of the "nucmer" and "promer" programs. Outputs
//            pairwise alignments to stdout. Works for both nucleotide and
//           amino-acid alignments.
//
//------------------------------------------------------------------------------

#include "delta.hh"
#include "tigrinc.hh"
#include "translate.hh"
#include "sw_alignscore.hh"
#include <vector>
#include <algorithm>
using namespace std;

//-- Output this many sequence characters per line
#define CHARS_PER_LINE 60

//------------------------------------------------------------- Constants ----//

const char NUCMER_MISMATCH_CHAR = '^';
const char NUCMER_MATCH_CHAR = ' ';
const char PROMER_SIM_CHAR = '+';
const char PROMER_MISMATCH_CHAR = ' ';

//-- Note: if coord exceeds LINE_PREFIX_LEN - 1 digits,
//         increase these accordingly
#define LINE_PREFIX_LEN 11
#define PREFIX_FORMAT "%-10ld "

#define DEFAULT_SCREEN_WIDTH 60
int Screen_Width = DEFAULT_SCREEN_WIDTH;



//------------------------------------------------------ Type Definitions ----//
struct AlignStats
     //-- Alignment statistics data structure
{
  long int sQ, eQ, sR, eR;              // start and end in Query and Reference
                                        // relative to the directional strand
  vector<long int> Delta;               // delta information
};



struct sR_Sort
//-- For sorting alignments by their sR coordinate
{
  bool operator( ) (const AlignStats & pA, const AlignStats & pB)
  {
    //-- sort sR
    if ( pA.sR < pB.sR )
      return true;
    else
      return false;
  }
};



struct sQ_Sort
//-- For sorting alignments by their sQ coordinate
{
  bool operator( ) (const AlignStats & pA, const AlignStats & pB)
  {
    //-- sort sQ
    if ( pA.sQ < pB.sQ )
      return true;
    else
      return false;
  }
};




//------------------------------------------------------ Global Variables ----//
bool isSortByQuery = false;              // -q option
bool isSortByReference = false;          // -r option



int DATA_TYPE = NUCMER_DATA;
int MATRIX_TYPE = BLOSUM62;

char InputFileName [MAX_LINE];
char RefFileName [MAX_LINE], QryFileName [MAX_LINE];

//------------------------------------------------- Function Declarations ----//
long int toFwd
     (long int coord, long int len, int frame);

void parseDelta
     (vector<AlignStats> & Aligns, char * IdR, char * IdQ);

void printAlignments
     (vector<AlignStats> Aligns, char * R, char * Q, char * IdR, char * IdQ);

void printHelp
     (const char * s);

void printUsage
     (const char * s);

long int revC
     (long int coord, long int len);



//-------------------------------------------------- Function Definitions ----//
int main
     (int argc, char ** argv)
{
  long int i;

  FILE * RefFile = NULL;
  FILE * QryFile = NULL;

  vector<AlignStats> Aligns;

  char * R, * Q;

  long int InitSize = INIT_SIZE;
  char Id [MAX_LINE], IdR [MAX_LINE], IdQ [MAX_LINE];

  //-- Parse the command line arguments
  {
    int ch, errflg = 0;
    optarg = NULL;

    while ( !errflg  &&  ((ch = getopt
                           (argc, argv, "hqro:c:w:x:")) != EOF) )
      switch (ch)
        {
        case 'h' :
	  printHelp (argv[0]);
	  exit (EXIT_SUCCESS);
          break;

	case 'q' :
	  isSortByQuery = true;
	  break;

	case 'r' :
	  isSortByReference = true;
	  break;

	case 'w' :
	  Screen_Width = atoi (optarg);
	  if ( Screen_Width <= LINE_PREFIX_LEN )
	    {
	      fprintf(stderr,
		      "WARNING: invalid screen width %d, using default\n",
		      DEFAULT_SCREEN_WIDTH);
	      Screen_Width = DEFAULT_SCREEN_WIDTH;
	    }
	  break;

	case 'x' :
	  MATRIX_TYPE = atoi (optarg);
	  if ( MATRIX_TYPE < 1 || MATRIX_TYPE > 3 )
	    {
	      fprintf(stderr,
		      "WARNING: invalid matrix type %d, using default\n",
		      MATRIX_TYPE);
	      MATRIX_TYPE = BLOSUM62;
	    }
	  break;

        default :
          errflg ++;
        }

    if ( errflg > 0  ||  argc - optind != 3 )
      {
        printUsage (argv[0]);
        exit (EXIT_FAILURE);
      }

    if ( isSortByQuery  &&  isSortByReference )
      fprintf (stderr,
               "WARNING: both -r and -q were passed, -q ignored\n");
  }

  strcpy (InputFileName, argv[optind ++]);
  strcpy (IdR, argv[optind ++]);
  strcpy (IdQ, argv[optind ++]);

  //-- Read in the alignment data
  parseDelta (Aligns, IdR, IdQ);

  //-- Find, and read in the reference sequence
  RefFile = File_Open (RefFileName, "r");
  InitSize = INIT_SIZE;
  R = (char *) Safe_malloc ( sizeof(char) * InitSize );
  while ( Read_String (RefFile, R, InitSize, Id, FALSE) )
    if ( strcmp (Id, IdR) == 0 )
      break;
  fclose (RefFile);
  if ( strcmp (Id, IdR) != 0 )
    {
      fprintf(stderr,"ERROR: Could not find %s in the reference file\n", IdR);
      exit (EXIT_FAILURE);
    }


  //-- Find, and read in the query sequence
  QryFile = File_Open (QryFileName, "r");
  InitSize = INIT_SIZE;
  Q = (char *) Safe_malloc ( sizeof(char) * InitSize );
  while ( Read_String (QryFile, Q, InitSize, Id, FALSE) )
    if ( strcmp (Id, IdQ) == 0 )
      break;
  fclose (QryFile);
  if ( strcmp (Id, IdQ) != 0 )
    {
      fprintf(stderr,"ERROR: Could not find %s in the query file\n", IdQ);
      exit (EXIT_FAILURE);
    }

  //-- Sort the alignment regions if user passed -r or -q option
  if ( isSortByReference )
    sort (Aligns.begin( ), Aligns.end( ), sR_Sort( ));
  else if ( isSortByQuery )
    sort (Aligns.begin( ), Aligns.end( ), sQ_Sort( ));


  //-- Output the alignments to stdout
  printAlignments (Aligns, R, Q,IdR,IdQ);

  return EXIT_SUCCESS;
}




long int toFwd
     (long int coord, long int len, int frame)

     // Switch relative coordinate to reference forward DNA strand

{
  long int newc = coord;

  if ( DATA_TYPE == PROMER_DATA )
    newc = newc * 3 - (3 - labs(frame));

  if ( frame < 0 )
    return revC ( newc, len );
  else
    return newc;
}




void parseDelta
     (vector<AlignStats> & Aligns, char * IdR, char * IdQ)

     // Read in the alignments from the desired region

{
  AlignStats aStats;                     //  single alignment region
  bool found = false;

  DeltaReader_t dr;
  dr.open (InputFileName);
  DATA_TYPE = dr.getDataType( ) == NUCMER_STRING ?
    NUCMER_DATA : PROMER_DATA;
  strcpy (RefFileName, dr.getReferencePath( ).c_str( ));
  strcpy (QryFileName, dr.getQueryPath( ).c_str( ));

  while ( dr.readNext( ) )
    {
      if ( dr.getRecord( ).idR == IdR  &&
	   dr.getRecord( ).idQ == IdQ )
	{
	  found = true;
	  break;
	}
    }
  if ( !found )
    {
      fprintf(stderr, "ERROR: Could not find any alignments for %s and %s\n",
	      IdR, IdQ);
      exit (EXIT_FAILURE);
    }

  for ( unsigned int i = 0; i < dr.getRecord( ).aligns.size( ); i ++ )
    {
      aStats.sR = dr.getRecord( ).aligns[i].sR;
      aStats.eR = dr.getRecord( ).aligns[i].eR;
      aStats.sQ = dr.getRecord( ).aligns[i].sQ;
      aStats.eQ = dr.getRecord( ).aligns[i].eQ;

      aStats.Delta = dr.getRecord( ).aligns[i].deltas;

      //-- Add the new alignment
      Aligns.push_back (aStats);
    }
  dr.close( );

  return;
}

void printCigarChar(int matches,int mismatches, int insertions, int deletions, int skips){
  if(matches){
    assert(mismatches==0);
    assert(insertions==0);
    assert(deletions==0);
    assert(skips==0);
    printf("%dM",matches);
  }
  if(mismatches){
    assert(matches==0);
    assert(insertions==0);
    assert(deletions==0);
    assert(skips==0);
    printf("%dX",mismatches);
  }
  if(insertions){
    assert(matches==0);
    assert(mismatches==0);
    assert(deletions==0);
    assert(skips==0);
    printf("%dI",insertions);
  }
  if(deletions){
    assert(matches==0);
    assert(mismatches==0);
    assert(insertions==0);
    assert(skips==0);
    printf("%dD",deletions);
  }
  if(skips){
    assert(matches==0);
    assert(mismatches==0);
    assert(insertions==0);
    assert(deletions==0);
    printf("%dS",skips);
  }
}		      


void printAlignments
     (vector<AlignStats> Aligns, char * R, char * Q, char * IdR, char * IdQ)

     // Print the alignments to the screen

{
  vector<AlignStats>::iterator Ap;
  vector<long int>::iterator Dp;
  int index = 1;
  
  char * A[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  char * B[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  int Ai, Bi, i;

  int Sign;
  long int Delta;
  long int Total, Errors, Remain;

  long int sR, eR, sQ, eQ;
  long int Apos, Bpos;
  long int SeqLenR, SeqLenQ;
  int frameR, frameQ;

  int matches=0;
  int mismatches=0;
  int skips=0;
  int insertions=0;
  int deletions=0;

  int ct = 0;

  //Set sequence lengths
  SeqLenR = strlen (R + 1);
  SeqLenQ = strlen (Q + 1);

  //Store sequence
  if ( DATA_TYPE == NUCMER_DATA )
    {
      A[1] = R;
      A[4] = (char *) Safe_malloc ( sizeof(char) * (SeqLenR + 2) );
      strcpy ( A[4] + 1, A[1] + 1 );
      A[4][0] = '\0';
      Reverse_Complement ( A[4], 1, SeqLenR );

      B[1] = Q;
      B[4] = (char *) Safe_malloc ( sizeof(char) * (SeqLenQ + 2) );
      strcpy ( B[4] + 1, B[1] + 1 );
      B[4][0] = '\0';
      Reverse_Complement ( B[4], 1, SeqLenQ );
    }

  for ( Ap = Aligns.begin( ); Ap < Aligns.end( ); Ap ++ )
    { 
      printf ("%s %s %d %d %d %d ",IdR,IdQ, Ap->sR,Ap->eR,Ap->sQ,Ap->eQ);
      index++;
      ct = 0;
      sR = Ap->sR;
      eR = Ap->eR;
      sQ = Ap->sQ;
      eQ = Ap->eQ;
      //-- Get the coords and frame right
      frameR = 1;
      if ( sR > eR )
	{
	  sR = revC (sR, SeqLenR);
	  eR = revC (eR, SeqLenR);
	  frameR += 3;
	}
      frameQ = 1;
      if ( sQ > eQ )
	{
	  sQ = revC (sQ, SeqLenQ);
	  eQ = revC (eQ, SeqLenQ);
	 
	  frameQ += 3;
	}

      Ai = frameR;
      Bi = frameQ;
      if ( frameR > 3 )
	frameR = -(frameR - 3);
      if ( frameQ > 3 )
	frameQ = -(frameQ - 3);

      //      skips = Ap->sR-1;

      Apos = sR;
      Bpos = sQ;

      Errors = 0;
      Total = 0;
      Remain = eR - sR + 1;

      for ( Dp = Ap->Delta.begin( );
	    Dp < Ap->Delta.end( ) &&
	    *Dp != 0; Dp ++ )
	{ 

	  Delta = *Dp;
	  Sign = Delta > 0 ? 1 : -1;
	  Delta = labs ( Delta );


	  //-- For all the bases before the next indel
	  for ( i = 1; i < Delta; i ++ )
	    {
	      if(A[Ai][Apos] == B[Bi][Bpos]){
		if(matches){
		  matches++;
		}
		else{
		  printCigarChar(matches,mismatches,insertions,deletions,skips);
		  matches=0;
		  mismatches=0;
		  insertions=0;
		  deletions=0;
		  skips=0;
		  matches++;
		}
		//		fprintf(Output,"%c",A[Ai][Apos]);
		if ( ++ ct == CHARS_PER_LINE ){
		  ct = 0;
		  //		  fprintf(Output, "\n");
		}
	      }
	      else{
		if(mismatches){
		  mismatches++;
		}
		else{
		  printCigarChar(matches,mismatches,insertions,deletions,skips);
		  matches=0;
		  mismatches=0;
		  insertions=0;
		  deletions=0;
		  skips=0;
		  mismatches++;
		}
		//fprintf(Output,"%c",A[Ai][Apos]);
		if ( ++ ct == CHARS_PER_LINE ){
		  ct = 0;
		  //fprintf(Output, "\n");
		}
	      }	
	      Apos ++;
	      Bpos ++;
	    }
	  //-- For the indel
	  Remain -= i - 1;
	  
	  if ( Sign == 1 ) {
	    if(insertions){
	      insertions++;
	    }
	    else{
	      printCigarChar(matches,mismatches,insertions,deletions,skips);
	      matches=0;
	      mismatches=0;
	      insertions=0;
	      deletions=0;
	      skips=0;
	      insertions++;
	    }
	    //fprintf(Output,"%c",A[Ai][Apos]);
	    if ( ++ ct == CHARS_PER_LINE ){
	      ct = 0;
	      //fprintf(Output, "\n");
	    }
	    Apos ++;
	    Remain --;
	  }
	  else {
	    if(deletions){
	      deletions++;
	    }
	    else{
	      printCigarChar(matches,mismatches,insertions,deletions,skips);
	      matches=0;
	      mismatches=0;
	      insertions=0;
	      deletions=0;
	      skips=0;
	      deletions++;
 	    }
	    Bpos ++;
	    Total ++;
	  }
	}
      //-- For all the bases remaining after the last indel
      for ( i = 0; i < Remain; i ++ )
	{
	  if(A[Ai][Apos] == B[Bi][Bpos]){
	    if(matches){
	      matches++;
	    }
	    else{
	      printCigarChar(matches,mismatches,insertions,deletions,skips);
	      matches=0;
	      mismatches=0;
	      insertions=0;
	      deletions=0;
	      skips=0;
	      matches++;
	    }
	    //fprintf(Output,"%c",A[Ai][Apos]);
	    if ( ++ ct == CHARS_PER_LINE ){
	      ct = 0;
	      //fprintf(Output, "\n");
	    }
	  }
	  else{
	    if(mismatches){
	      mismatches++;
	    }
	    else{
	      printCigarChar(matches,mismatches,insertions,deletions,skips);
	      matches=0;
	      mismatches=0;
	      insertions=0;
	      deletions=0;
	      skips=0;
	      mismatches++;
	    }
	    //fprintf(Output,"%c",A[Ai][Apos]);
	    if ( ++ ct == CHARS_PER_LINE ){
	      ct = 0;
	      //fprintf(Output, "\n");
	    }
	  }
	  Apos ++;
	  Bpos ++;
	}
      printCigarChar(matches,mismatches,insertions,deletions,skips);
      matches=0;
      mismatches=0;
      insertions=0;
      deletions=0;
      skips=0;
      //fprintf(Output, "\n");
      printf ("\n");
    }
}
  

void printHelp
     (const char * s)

      //  Display the program's help information to stderr

{
  fprintf (stderr,
           "\nUSAGE: %s  [options]  <deltafile>  <ref ID>  <qry ID>\n\n", s);
  fprintf (stderr,
       "-h            Display help information\n"
       "-q            Sort alignments by the query start coordinate\n"
       "-r            Sort alignments by the reference start coordinate\n"
       "-w int        Set the screen width - default is 60\n"
       "-x int        Set the matrix type - default is 2 (BLOSUM 62),\n"
       "              other options include 1 (BLOSUM 45) and 3 (BLOSUM 80)\n"
       "              note: only has effect on amino acid alignments\n\n");
  fprintf (stderr,
       "  Input is the .delta output of either the \"nucmer\" or the\n"
       "\"promer\" program passed on the command line.\n"
       "  Output is to stdout, and consists of all the alignments between the\n"
       "query and reference sequences identified on the command line.\n"
       "  NOTE: No sorting is done by default, therefore the alignments\n"
       "will be ordered as found in the <deltafile> input.\n\n");
  return;
}




void printUsage
     (const char * s)

      //  Display the program's usage information to stderr.

{
  fprintf (stderr,
           "\nUSAGE: %s  [options]  <deltafile>  <ref ID>  <qry ID>\n\n", s);
  fprintf (stderr, "Try '%s -h' for more information.\n", s);
  return;
}




long int revC
     (long int coord, long int len)
{
  return len - coord + 1;
}
