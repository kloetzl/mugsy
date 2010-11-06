#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/seeds.h>

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

//template <typename TAlgorithmSpec>
void Test_SimpleSeeds()
{

//____________________________________________________________________________
// Standard Functions
	Seed<int,SimpleSeed> seed;
	Seed<int, SimpleSeed> seed1(0,0,7);
	Seed<int, SimpleSeed> seed2(0,1,7,5);

	SEQAN_TASSERT(startDiagonal(seed1)==0);
	SEQAN_TASSERT(endDiagonal(seed1)==0);
	SEQAN_TASSERT(startDiagonal(seed2)==1);
	SEQAN_TASSERT(endDiagonal(seed2)==-2);
	SEQAN_TASSERT(leftDim0(seed2)==0);
	SEQAN_TASSERT(rightDim0(seed2)==7);
	SEQAN_TASSERT(leftDim1(seed2)==1);
	SEQAN_TASSERT(rightDim1(seed2)==5);
	SEQAN_TASSERT(leftDiagonal(seed2)==1);
	SEQAN_TASSERT(rightDiagonal(seed2)==-2);
	SEQAN_TASSERT(length(seed2) == 8);

	SEQAN_TASSERT(leftPosition(seed2,0) == 0);
	SEQAN_TASSERT(leftPosition(seed2,1) == 1);
	SEQAN_TASSERT(rightPosition(seed2,0) == 7);
	SEQAN_TASSERT(rightPosition(seed2,1) == 5);

	setLeftDim0(seed2,3);
	setRightDim0(seed2,9);
	setLeftDim1(seed2,5);
	setRightDim1(seed2,12);
	setLeftDiagonal(seed2,29);
	setRightDiagonal(seed2,7);
	SEQAN_TASSERT(leftDim0(seed2)==3);
	SEQAN_TASSERT(rightDim0(seed2)==9);
	SEQAN_TASSERT(leftDim1(seed2)==5);
	SEQAN_TASSERT(rightDim1(seed2)==12);
	SEQAN_TASSERT(leftDiagonal(seed2)==29);
	SEQAN_TASSERT(rightDiagonal(seed2)==7);

//____________________________________________________________________________
// Merge Algorithms
	Seed<int, SimpleSeed> seed3(0,0,7);
	Seed<int, SimpleSeed> seed4(4,5,7);
	_mergeTwoSeeds(seed3,seed4,Merge());
	SEQAN_TASSERT(startDiagonal(seed3)==0);
	SEQAN_TASSERT(endDiagonal(seed3)==1);
	SEQAN_TASSERT(startDiagonal(seed3)==0);
	SEQAN_TASSERT(endDiagonal(seed3)==1);
	SEQAN_TASSERT(leftDim0(seed3)==0);
	SEQAN_TASSERT(rightDim0(seed3)==10);
	SEQAN_TASSERT(leftDim1(seed3)==0);
	SEQAN_TASSERT(rightDim1(seed3)==11);
	SEQAN_TASSERT(leftDiagonal(seed3)==1);
	SEQAN_TASSERT(rightDiagonal(seed3)==0);
	SEQAN_TASSERT(length(seed3) == 11);

	
	Seed<int, SimpleSeed> seed5(0,0,7);
	_mergeTwoSeeds(seed5,4,5,7,Merge());
	SEQAN_TASSERT(startDiagonal(seed5)==0);
	SEQAN_TASSERT(endDiagonal(seed5)==1);
	SEQAN_TASSERT(startDiagonal(seed5)==0);
	SEQAN_TASSERT(endDiagonal(seed5)==1);
	SEQAN_TASSERT(leftDim0(seed5)==0);
	SEQAN_TASSERT(rightDim0(seed5)==10);
	SEQAN_TASSERT(leftDim1(seed5)==0);
	SEQAN_TASSERT(rightDim1(seed5)==11);
	SEQAN_TASSERT(leftDiagonal(seed5)==1);
	SEQAN_TASSERT(rightDiagonal(seed5)==0);
	SEQAN_TASSERT(length(seed5) == 11);

	Seed<int, SimpleSeed> seed6(0,0,7);
	_mergeTwoSeeds(seed6,4,5,10,11,Merge());
	SEQAN_TASSERT(startDiagonal(seed6)==0);
	SEQAN_TASSERT(endDiagonal(seed6)==1);
	SEQAN_TASSERT(startDiagonal(seed6)==0);
	SEQAN_TASSERT(endDiagonal(seed6)==1);
	SEQAN_TASSERT(leftDim0(seed6)==0);
	SEQAN_TASSERT(rightDim0(seed6)==10);
	SEQAN_TASSERT(leftDim1(seed6)==0);
	SEQAN_TASSERT(rightDim1(seed6)==11);
	SEQAN_TASSERT(leftDiagonal(seed6)==1);
	SEQAN_TASSERT(rightDiagonal(seed6)==0);
	SEQAN_TASSERT(length(seed6) == 11);

	Score<int,Simple> matrix(2,-1,-1);
	Seed<int, SimpleSeed> seed10(0,0,7);
	int score10 = 14;
	Seed<int, SimpleSeed> seed11(0,0,7);
	int score11 = 14;
	Seed<int, SimpleSeed> seed12(0,0,7);
	int score12 = 14;
	Seed<int, SimpleSeed> seed13(4,5,7);
	int score13 = 14;


	_mergeTwoSeedsScore(seed10, score10, seed13, score13, matrix, Manhattan(), Merge());
	SEQAN_TASSERT(leftDim0(seed10)==0);
	SEQAN_TASSERT(rightDim0(seed10)==10);
	SEQAN_TASSERT(score10 = 20);
	_mergeTwoSeedsScore(seed11, score11, 4, 5, 7, score13, matrix, Manhattan(), Merge());
	SEQAN_TASSERT(leftDim0(seed11)==0);
	SEQAN_TASSERT(rightDim0(seed11)==10);
	SEQAN_TASSERT(score11 = 20);         
	_mergeTwoSeedsScore(seed12, score12, 4, 5, 10, 11, score13, matrix, Manhattan(), Merge());
	SEQAN_TASSERT(leftDim0(seed12)==0);
	SEQAN_TASSERT(rightDim0(seed12)==10);
	SEQAN_TASSERT(score12 = 20);

//____________________________________________________________________________
// Extension Algorithms



	String<Dna> query =	   "AAACCCTTTGGGTTTTT";
	String<Dna> database = "AACCCCTTTGGTGAAAAA";
	Seed<int, SimpleSeed> seed7(4,4,3);
	extendSeed(seed7, query, database, 2, MatchExtend());
	SEQAN_TASSERT(leftDim0(seed7)==3);
	SEQAN_TASSERT(rightDim0(seed7)==10);
	SEQAN_TASSERT(leftDim1(seed7)==3);
	SEQAN_TASSERT(rightDim1(seed7)==10);

	Seed<int, SimpleSeed> seed8(4,4,3);
	extendSeed(seed8, 2, matrix, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(leftDim0(seed8)==0);
	SEQAN_TASSERT(rightDim0(seed8)==10);
	SEQAN_TASSERT(leftDim1(seed8)==0);
	SEQAN_TASSERT(rightDim1(seed8)==10);

	Seed<int, SimpleSeed> seed9(4,4,3);
	extendSeed(seed9, 1, matrix, query, database, 2, GappedXDrop());
	SEQAN_TASSERT(leftDim0(seed9)==0);
	SEQAN_TASSERT(rightDim0(seed9)==12);
	SEQAN_TASSERT(leftDim1(seed9)==0);
	SEQAN_TASSERT(rightDim1(seed9)==13);
}

void Test_MultiSeeds(){

//____________________________________________________________________________
// Standard Functions
	Seed<int,ChainedSeed> seed;
	Seed<int, ChainedSeed> seed1(4,5,7);
	SEQAN_TASSERT(startDiagonal(seed1)==1);
	SEQAN_TASSERT(endDiagonal(seed1)==1);
	SEQAN_TASSERT(leftDim0(seed1)==4);
	SEQAN_TASSERT(rightDim0(seed1)==10);
	SEQAN_TASSERT(leftDim1(seed1)==5);
	SEQAN_TASSERT(rightDim1(seed1)==11);
	SEQAN_TASSERT(leftDiagonal(seed1)==1);
	SEQAN_TASSERT(rightDiagonal(seed1)==1);
	SEQAN_TASSERT(length(seed1) == 7);
	SEQAN_TASSERT(_getFirstDiag(seed1).i2 == 5);

	Seed<int, ChainedSeed> const seed132(4,5,7);
	SEQAN_TASSERT(_getLastDiag(seed132).i1 == 4);


	
	Seed<int,ChainedSeed> seed2(0,0,0);
	setLeftDim0(seed2,3);
	setRightDim0(seed2,9);
	setLeftDim1(seed2,2);
	setRightDim1(seed2,12);
	setLeftDiagonal(seed2,29);
	setRightDiagonal(seed2,7);
	SEQAN_TASSERT(leftDim0(seed2)==2);
	SEQAN_TASSERT(rightDim0(seed2)==12);
	SEQAN_TASSERT(leftDim1(seed2)==2);
	SEQAN_TASSERT(rightDim1(seed2)==12);
	SEQAN_TASSERT(leftDiagonal(seed2)==29);
	SEQAN_TASSERT(rightDiagonal(seed2)==7);

//____________________________________________________________________________
// Merge Algorithms
	Seed<int, ChainedSeed> seed3(0,0,7);
	Seed<int, ChainedSeed> seed4(4,5,7);
	_mergeTwoSeeds(seed3,seed4,Merge());
	SEQAN_TASSERT(startDiagonal(seed3)==0);
	SEQAN_TASSERT(endDiagonal(seed3)==1);
	SEQAN_TASSERT(startDiagonal(seed3)==0);
	SEQAN_TASSERT(endDiagonal(seed3)==1);
	SEQAN_TASSERT(leftDim0(seed3)==0);
	SEQAN_TASSERT(rightDim0(seed3)==10);
	SEQAN_TASSERT(leftDim1(seed3)==0);
	SEQAN_TASSERT(rightDim1(seed3)==11);
	SEQAN_TASSERT(leftDiagonal(seed3)==1);
	SEQAN_TASSERT(rightDiagonal(seed3)==0);
	SEQAN_TASSERT(length(seed3) == 11);

	
	Seed<int, ChainedSeed> seed5(0,0,7);
	_mergeTwoSeeds(seed5,4,5,7,Merge());
	SEQAN_TASSERT(startDiagonal(seed5)==0);
	SEQAN_TASSERT(endDiagonal(seed5)==1);
	SEQAN_TASSERT(startDiagonal(seed5)==0);
	SEQAN_TASSERT(endDiagonal(seed5)==1);
	SEQAN_TASSERT(leftDim0(seed5)==0);
	SEQAN_TASSERT(rightDim0(seed5)==10);
	SEQAN_TASSERT(leftDim1(seed5)==0);
	SEQAN_TASSERT(rightDim1(seed5)==11);
	SEQAN_TASSERT(leftDiagonal(seed5)==1);
	SEQAN_TASSERT(rightDiagonal(seed5)==0);
	SEQAN_TASSERT(length(seed5) == 11);


	Score<int,Simple> matrix(2,-1,-1);

    Seed<int, ChainedSeed> seed10(0,0,7);
	int score10 = 14;
	Seed<int, ChainedSeed> seed11(0,0,7);
	int score11 = 14;
	Seed<int, ChainedSeed> seed12(4,5,7);
	int score12 = 14;


	_mergeTwoSeedsScore(seed10, score10, 4, 5, 7, 14, matrix, Manhattan(), Merge());
	_mergeTwoSeedsScore(seed10, score10, 3, 4, 10, 20, matrix, Manhattan(), Merge());
    SEQAN_TASSERT(leftDim0(seed10)==0);
	SEQAN_TASSERT(rightDim0(seed10)==12);
	SEQAN_TASSERT(score10 = 26);
	
	_mergeTwoSeedsScore(seed11, score11, seed12, score12, matrix, Manhattan(), Merge());
	SEQAN_TASSERT(leftDim0(seed11)==0);
	SEQAN_TASSERT(rightDim0(seed11)==10);
	SEQAN_TASSERT(score11 = 20);

//____________________________________________________________________________
// Extension Algorithms

	String<Dna> query =	   "AAACCCTTTGGGTTTTT";
	String<Dna> database = "AACCCCTTTGGTGAAAAA";
	Seed<int, ChainedSeed> seed7(4,4,3);
	extendSeed(seed7, query, database, 2, MatchExtend());
	SEQAN_TASSERT(leftDim0(seed7)==3);
	SEQAN_TASSERT(rightDim0(seed7)==10);
	SEQAN_TASSERT(leftDim1(seed7)==3);
	SEQAN_TASSERT(rightDim1(seed7)==10);

	Seed<int, ChainedSeed> seed8(4,4,3);
	extendSeed(seed8, 2, matrix, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(leftDim0(seed8)==0);
	SEQAN_TASSERT(rightDim0(seed8)==10);
	SEQAN_TASSERT(leftDim1(seed8)==0);
	SEQAN_TASSERT(rightDim1(seed8)==10);

	Seed<int, ChainedSeed> seed9(4,4,3);
	extendSeed(seed9, 1, matrix, query, database, 2, GappedXDrop());
	SEQAN_TASSERT(leftDim0(seed9)==0);
	SEQAN_TASSERT(rightDim0(seed9)==12);
	SEQAN_TASSERT(leftDim1(seed9)==0);
	SEQAN_TASSERT(rightDim1(seed9)==13);

//____________________________________________________________________________
// Alignment Calculation
	Align<String<Dna>, ArrayGaps> aligned;
	getAlignment(seed8, aligned, query, database, matrix);
	SEQAN_TASSERT(row(aligned,0)== "AAACCCTTTGG");
	SEQAN_TASSERT(row(aligned,1)== "AACCCCTTTGG");
	
	
	
//____________________________________________________________________________
// Score Calculation
	SEQAN_TASSERT(scoreSeed(seed8, query, database, matrix) == 19);
}

void Main_Seeds(){
	SEQAN_TREPORT("TEST BEGIN")
	Test_SimpleSeeds();
	Test_MultiSeeds();
	debug::verifyCheckpoints("projects/library/seqan/seeds/seed_base.h");
	debug::verifyCheckpoints("projects/library/seqan/seeds/seed_multi.h");
	SEQAN_TREPORT("TEST END")
}
