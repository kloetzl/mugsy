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
#include <map>
#include <string>
using namespace std;

//------------------------------------------------------------- Constants ----//

const char NUCMER_MISMATCH_CHAR = '^';
const char NUCMER_MATCH_CHAR = ' ';
const char PROMER_SIM_CHAR = '+';
const char PROMER_MISMATCH_CHAR = ' ';

//-- Note: if coord exceeds LINE_PREFIX_LEN - 1 digits,
//         increase these accordingly
#define LINE_PREFIX_LEN 11
#define PREFIX_FORMAT "%-10ld "

#define DEFAULT_SCREEN_WIDTH 100000
int Screen_Width = DEFAULT_SCREEN_WIDTH;



//------------------------------------------------------ Type Definitions ----//
struct AlignStats
     //-- Alignment statistics data structure
{
  long int sQ, eQ, sR, eR;              // start and end in Query and Reference
                                        // relative to the directional strand
  vector<long int> Delta;               // delta information
  std::string idR;         //!< reference contig ID
  std::string idQ;         //!< query contig ID
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
bool forceDNA = true;          // show DNA alignments even for promer generated alignments

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
(vector<AlignStats> Aligns, char * R, char * Q,  map<string,char *> & seqsMap);

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
  IdR[0]='\0';
  IdQ[0]='\0';

  map<string, char *> seqsMap;

  //-- Parse the command line arguments
  {
    int ch, errflg = 0;
    optarg = NULL;

    while ( !errflg  &&  ((ch = getopt
                           (argc, argv, "hqrw:x:")) != EOF) )
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

    if ( errflg > 0)
      {
        printUsage (argv[0]);
        exit (EXIT_FAILURE);
      }

    if ( isSortByQuery  &&  isSortByReference )
      fprintf (stderr,
               "WARNING: both -r and -q were passed, -q ignored\n");
  }

  strcpy (InputFileName, argv[optind ++]);
  if((argc - optind) >=1)
    strcpy (IdR, argv[optind ++]);
  if((argc - optind) >=1)
    strcpy (IdQ, argv[optind ++]);

  //-- Read in the alignment data
  parseDelta (Aligns, IdR, IdQ);

  //-- Find, and read in the reference sequence
  RefFile = File_Open (RefFileName, "r");
  InitSize = INIT_SIZE;
  Read_File (RefFile, InitSize, seqsMap, FALSE);
  //printf("Seqmap size %d\n",seqsMap.size());
  fclose (RefFile);
  //-- Find, and read in the query sequence
  QryFile = File_Open (QryFileName, "r");
  InitSize = INIT_SIZE;
  Read_File (QryFile, InitSize, seqsMap, FALSE);
  //printf("Seqmap size %d\n",seqsMap.size());
  fclose (QryFile);

  //-- Sort the alignment regions if user passed -r or -q option
  if ( isSortByReference )
    sort (Aligns.begin( ), Aligns.end( ), sR_Sort( ));
  else if ( isSortByQuery )
    sort (Aligns.begin( ), Aligns.end( ), sQ_Sort( ));


  //-- Output the alignments to stdout
  //  printf("%s %s\n\n", RefFileName, QryFileName);
  //for ( i = 0; i < Screen_Width; i ++ ) printf("=");
  //printf("\n-- Alignments between %s and %s\n\n", IdR, IdQ);s
  printf("##maf version=1 scoring=single_cov2\n");
  printAlignments (Aligns, R, Q, seqsMap);
  //  printf("\n");
  //for ( i = 0; i < Screen_Width; i ++ ) printf("=");
  //printf("\n");

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
  bool foundany = false;

  DeltaReader_t dr;
  dr.open (InputFileName);

  if(forceDNA)
    DATA_TYPE = NUCMER_DATA;
  else
    DATA_TYPE = dr.getDataType( ) == NUCMER_STRING ?
      NUCMER_DATA : PROMER_DATA;

  strcpy (RefFileName, dr.getReferencePath( ).c_str( ));
  strcpy (QryFileName, dr.getQueryPath( ).c_str( ));

  //  while ( dr.readNext( ) )
  //{
  //  if(IdR != NULL && IdQ != NULL 
  // && dr.getRecord( ).idR == IdR  &&
  // dr.getRecord( ).idQ == IdQ )
  //{
  //  found = true;
  //  break;
  //}
  //  else 
  //{ 
  //  if(IdR != NULL && dr.getRecord( ).idR == IdR ){
  //    found = true;
  //    break;
  //  }
  //  else
  //    {
  //      if(IdQ != NULL && dr.getRecord( ).idQ == IdQ ){
  //	found = true;
  //	break;
  //      }
  //    }
  //}
  //}

  //printf ("IdR:%s %d IdQ:%s %d\n",IdR,strlen(IdR),IdQ,strlen(IdQ));
  while ( dr.readNext( ) ){
    for ( unsigned int i = 0; i < dr.getRecord( ).aligns.size( ); i ++ )
      {
	if(strlen(IdR) == 0 && strlen(IdQ) == 0){
	  found=true;
	}
	else{
	  if(strlen(IdR) != 0 && strlen(IdQ) != 0){
	    if(dr.getRecord( ).idR == IdR && dr.getRecord( ).idQ == IdQ){
	      found=true;
	    }
	  }
	  else{
	    if(strlen(IdR) != 0 && dr.getRecord( ).idR == IdR){
	      found=true;
	    }
	    else{
	      if(strlen(IdQ) != 0 && dr.getRecord( ).idQ == IdQ){
		found=true;
	      }
	    }
	  }
	}
	if(found){
	  aStats.sR = dr.getRecord( ).aligns[i].sR;
	  aStats.eR = dr.getRecord( ).aligns[i].eR;
	  aStats.sQ = dr.getRecord( ).aligns[i].sQ;
	  aStats.eQ = dr.getRecord( ).aligns[i].eQ;
	  aStats.idR = dr.getRecord( ).idR;
	  aStats.idQ = dr.getRecord( ).idQ;
	  //printf("Saving match ref=%s query=%s\n",aStats.idR.c_str(),aStats.idQ.c_str());
	  aStats.Delta = dr.getRecord( ).aligns[i].deltas;
	
	  //-- Add the new alignment
	  Aligns.push_back (aStats);
	  foundany=true;
	}
      }
    found=false;
  }

  dr.close( );

  if ( !foundany )
   {
     fprintf(stderr, "ERROR: Could not find any alignments for %s and %s\n",
  	      IdR, IdQ);
     printf("##maf version=1 scoring=single_cov2\n");
     
     printf("##eof maf");
      
      exit (EXIT_FAILURE);
    }

  return;
}




void printAlignments
(vector<AlignStats> Aligns, char * R, char * Q, map<string, char *> & seqsMap)

     // Print the alignments to the screen

{

  const char * IdR;
  const char * IdQ;

  map<string, char *>::iterator finditer;

  map<pair<string,int>, char *> seqsMapArray;
  map<pair<string,int>, char *>::iterator seqsiter;

  vector<AlignStats>::iterator Ap;
  vector<long int>::iterator Dp;

  char * A[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  char * B[7] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL};
  int Ai, Bi, i;

  char Buff1 [Screen_Width + 1],
    Buff2 [Screen_Width + 1];
  //Buff3 [Screen_Width + 1];

  int Sign;
  long int Delta;
  long int Total, Errors, Remain;
  long int Pos;

  long int sR, eR, sQ, eQ;
  long int Apos, Bpos;
  long int SeqLenR, SeqLenQ;
  int frameR, frameQ;

  //for ( i = 0; i < LINE_PREFIX_LEN; i ++ )
  //Buff3[i] = ' ';
  for ( Ap = Aligns.begin( ); Ap < Aligns.end( ); Ap ++ )
    {	
      //HACK, shortcut to test perf
      //memset(&Buff1,'Z',Screen_Width + 1);
      //memset(&Buff2,'Z',Screen_Width + 2);

      sR = Ap->sR;
      eR = Ap->eR;
      sQ = Ap->sQ;
      eQ = Ap->eQ;
      IdR = Ap->idR.c_str();
      IdQ = Ap->idQ.c_str();

      finditer = seqsMap.find(Ap->idR);
      //printf("Looking for R:\"%s\" in map of size %d\n",IdR,seqsMap.size());
      assert(finditer != seqsMap.end());
      R = finditer->second;
      SeqLenR = strlen(R+1);

      if(DATA_TYPE == NUCMER_DATA){
	seqsiter = seqsMapArray.find(make_pair(Ap->idR,1));
	if(seqsiter == seqsMapArray.end()){
	  A[1] = R;
	  A[4] = (char *) Safe_malloc ( sizeof(char) * (SeqLenR + 2) );
	  strcpy ( A[4] + 1, A[1] + 1 );
	  A[4][0] = '\0';
	  Reverse_Complement ( A[4], 1, SeqLenR );
	  //printf("#Allocating memory for %s\n",Ap->idR.c_str());
	  seqsMapArray.insert(make_pair(make_pair(Ap->idR,1),A[1]));
	  seqsMapArray.insert(make_pair(make_pair(Ap->idR,4),A[4]));
	}
	else{
	  A[1] = seqsiter->second;
	  seqsiter = seqsMapArray.find(make_pair(Ap->idR,4));
	  assert(seqsiter != seqsMapArray.end());
	  A[4] = seqsiter->second;
	}
      }
    
      finditer = seqsMap.find(Ap->idQ);
      //printf("Looking for Q:\"%s\" in map of size %d\n",IdQ,seqsMap.size());
      assert(finditer != seqsMap.end());
      Q = finditer->second;
      SeqLenQ = strlen(Q+1);

      if(DATA_TYPE == NUCMER_DATA){
	seqsiter = seqsMapArray.find(make_pair(Ap->idQ,1));
	if(seqsiter == seqsMapArray.end()){
	  B[1] = Q;
	  B[4] = (char *) Safe_malloc ( sizeof(char) * (SeqLenQ + 2) );
	  strcpy ( B[4] + 1, B[1] + 1 );
	  B[4][0] = '\0';
	  Reverse_Complement ( B[4], 1, SeqLenQ );
	  //printf("#Allocating memory for %s\n",Ap->idQ.c_str());
	  seqsMapArray.insert(make_pair(make_pair(Ap->idQ,1),B[1]));
	  seqsMapArray.insert(make_pair(make_pair(Ap->idQ,4),B[4]));
	}
	else{
	  B[1] = seqsiter->second;
	  //printf("#Looking for Q:\"%s\" in map of size %d\n",IdQ,seqsMapArray.size());
	  seqsiter = seqsMapArray.find(make_pair(Ap->idQ,4));
	  assert(seqsiter != seqsMapArray.end());
	  B[4] = seqsiter->second;
	}
      }

    
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
      
      if ( DATA_TYPE == PROMER_DATA )
	{
	  frameR += (sR + 2) % 3;
	  frameQ += (sQ + 2) % 3;
	  
	  //-- Translate the coordinates from DNA to Amino Acid
	  //   remeber that eR and eQ point to the last nucleotide in the codon
	  sR = (sR + 2) / 3;
	  eR = eR / 3;
	  sQ = (sQ + 2) / 3;
	  eQ = eQ / 3;
	}
      Ai = frameR;
      Bi = frameQ;
      if ( frameR > 3 )
	frameR = -(frameR - 3);
      if ( frameQ > 3 )
	frameQ = -(frameQ - 3);

      /*
      if ( A[Ai] == NULL ){
	assert ( DATA_TYPE == PROMER_DATA );
	A[Ai] = (char *) Safe_malloc ( sizeof(char) * ( SeqLenR / 3 + 2 ) );
	A[Ai][0] = '\0';
	Translate_DNA ( R, A[Ai], Ai );
      }

      if ( B[Bi] == NULL ){
	assert ( DATA_TYPE == PROMER_DATA );
	B[Bi] = (char *) Safe_malloc ( sizeof(char) * ( SeqLenQ / 3 + 2 ) );
	B[Bi][0] = '\0';
	Translate_DNA ( Q, B[Bi], Bi );
      }
      */
      //-- Generate the alignment
      printf("a score=%d\n",abs(Ap->eR - Ap->sR));


      //Loop over query and reference
      int query;
      for(query=0;query<2;query++){
	Apos = sR;
	Bpos = sQ;
	
	Errors = 0;
	Total = 0;
	Remain = eR - sR + 1;
	
	//sprintf(Buff1, PREFIX_FORMAT, toFwd (Apos, SeqLenR, frameR));
	//sprintf(Buff2, PREFIX_FORMAT, toFwd (Bpos, SeqLenQ, frameQ));
	Pos = 0;
	/*
	int rgaps=0;
	int qgaps=0;
	
	for ( Dp = Ap->Delta.begin( );
	      Dp < Ap->Delta.end( ) &&
		*Dp != 0; Dp ++ )
	  {
	    Delta = *Dp;
	    Sign = Delta > 0 ? 1 : -1;
	    Delta = labs ( Delta );
	    if(Sign < 0){
	      rgaps++;
	    }
	    else{
	      qgaps++;
	    }
	  }
	*/
	//printf("# gaps %d %d\n",rgaps,qgaps);
	if(query == 1){
	  //printf("#%s s:%d e:%d len:%d f:%d\n",IdQ,sQ,eQ,eQ-sQ+1,frameQ);
	  if(frameQ < 0){
	    //	  printf("s %s %d %d %s %d ", IdQ, SeqLenQ - (sQ-1+eQ-sQ+1), eQ-sQ+1, "-",SeqLenQ); 	
	    printf("s %s %d %d %s %d ", IdQ, sQ-1, eQ-sQ+1, "-",SeqLenQ); 	
	  }
	  else{
	    printf("s %s %d %d %s %d ", IdQ, sQ-1, eQ-sQ+1, "+",SeqLenQ); 	
	  }
	}
	else{ 
	  //printf("#%s s:%d e:%d len:%d f:%d\n",IdR,sR,eR,eR-sR+1,frameR);
	  if(frameR < 0){
	    //printf("s %s %d %d %s %d ", IdR, SeqLenR - (sR-1+eR-sR+1), eR-sR+1, "-",SeqLenR); 
	    printf("s %s %d %d %s %d ", IdR, sR-1, eR-sR+1, "-",SeqLenR); 
	  }
	  else{
	    printf("s %s %d %d %s %d ", IdR, sR-1, eR-sR+1, "+",SeqLenR); 
	  }
	}
	
	for ( Dp = Ap->Delta.begin( );
	      Dp < Ap->Delta.end( ) &&
		*Dp != 0; Dp ++ )
	  {
	    Delta = *Dp;
	    Sign = Delta > 0 ? 1 : -1;
	    Delta = labs ( Delta );
	    if(Pos+Delta-1 < Screen_Width){
	      if(query==0){
		memcpy(&Buff1[Pos],&A[Ai][Apos],Delta-1);
		//memset(&Buff1[Pos],'Z',Delta-1);
		Apos = Apos + Delta - 1;
	      }
	      else{
		memcpy(&Buff2[Pos],&B[Bi][Bpos],Delta-1);
		//memset(&Buff2[Pos],'Z',Delta-1);
		Bpos = Bpos + Delta - 1;
	      }
	      Pos = Pos + Delta - 1;
	      i = Delta;
	    }
	    else{
	      //-- For all the bases before the next indel
	      for ( i = 1; i < Delta; i ++ )
		{
		  if ( Pos >= Screen_Width )
		    {
		      if(query == 1){
			Buff2[Pos] = '\0';
			printf("%s", &Buff2);
		      }
		      else{
			Buff1[Pos] = '\0';
			printf("%s", &Buff1);
		      }
		      Pos = 0;
		    }
		  if(query==0){
		    Buff1[Pos] = A[Ai][Apos ++];
		  }
		  else{
		    Buff2[Pos] = B[Bi][Bpos ++];
		  }
		  Pos++;
		}
	    }

	    
	    //-- For the indel
	    Remain -= i - 1;
	    
	    if ( Pos >= Screen_Width )
	      {
		if(query == 1){
		  Buff2[Pos] = '\0';
		  printf("%s", &Buff2);
		}
		else{
		  Buff1[Pos] = '\0';
		  printf("%s", &Buff1);
		}
		Pos = 0;
	      }
		    
	    if ( Sign == 1 )
	      {
		if(query==0)
		  Buff1[Pos] = A[Ai][Apos ++];
		else
		  Buff2[Pos] = '-';
		Pos++;
		Remain --;
	      }
	    else
	      {
		if(query==0)
		  Buff1[Pos] = '-';
		else
		  Buff2[Pos] = B[Bi][Bpos ++];
		Pos++;
		Total ++;
	      }
	  }


	//-- For all the bases remaining after the last indel
	if(Pos+Remain < Screen_Width){
	  if(query==0){
	    memcpy(&Buff1[Pos],&A[Ai][Apos],Remain);
	    //memset(&Buff1[Pos],'Z',Remain);
	    Apos = Apos + Remain;
	  }
	  else{
	    memcpy(&Buff2[Pos],&B[Bi][Bpos],Remain);
	    //memset(&Buff2[Pos],'Z',Remain);
	    Bpos = Bpos + Remain;
	  }
	  Pos = Pos + Remain;
	}
	else{
	  for ( i = 0; i < Remain; i ++ )
	    {
	      if ( Pos >= Screen_Width )
		{
		  if(query == 1){
		    Buff2[Pos] = '\0';
		    printf("%s", &Buff2);
		  }
		  else{
		    Buff1[Pos] = '\0';
		    printf("%s", &Buff1);
		  }
		  Pos = 0;
		}
	      if(query==0)
		Buff1[Pos] = A[Ai][Apos ++];
	      else
		Buff2[Pos] = B[Bi][Bpos ++];
	      Pos++;
	    }
	}

	
	//-- For the remaining buffered 
	if ( Pos > 0)
	  {
	    if(query == 1){
	      Buff2[Pos] = '\0';
	      printf("%s", &Buff2);
	    }
	    else{
	      Buff1[Pos] = '\0';
	      printf("%s", &Buff1);
	    }
	    Pos = 0;
	  }
	printf("\n");
	if(query==1){
	  printf("\n");
	}
      }
    }
  //SVA, leaks here because I'm saving the all the seqs
  //-- Free the sequences, except for the originals
  for ( i = 0; i < 7; i ++ )
    {
      if ( (DATA_TYPE != NUCMER_DATA || i != 1)  &&  A[i] != NULL )
	free ( A[i] );
      if ( (DATA_TYPE != NUCMER_DATA || i != 1)  &&  B[i] != NULL )
	free ( B[i] );
    }
  printf("##eof maf");
  return;
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
