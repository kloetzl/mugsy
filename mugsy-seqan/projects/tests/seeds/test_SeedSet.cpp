#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/seeds.h>


using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

//template <typename TAlgorithmSpec>
void Test_SeedSet_base()
{
	SeedSet<int,SimpleSeed, DefaultNoScore> set1;
	SeedSet<int,ChainedSeed, DefaultNoScore> set2;
	SeedSet<int,SimpleSeed, DefaultNoScore> set3(100,0);
	SeedSet<int,ChainedSeed, DefaultNoScore> set4(100,0);

	SEQAN_TASSERT(qualityValue(set3)==0);
	SEQAN_TASSERT(maximumDistance(set3)==100);
	setQualityValue(set3,1);
	setMaximumDistance(set3,99);
	SEQAN_TASSERT(qualityValue(set3)==1);
	SEQAN_TASSERT(maximumDistance(set3)==99);
	SEQAN_TASSERT(length(set3)==0);	
	//SEQAN_TASSERT(capacity(set3)==1);

	Seed<int,SimpleSeed> seed(0,0,3);
	Seed<int,SimpleSeed> seed2(3,4,3);
	appendValue(set3,seed);
	appendValue(set3,seed2);
	SEQAN_TASSERT(length(set3)==2);	

	SeedSet<int,SimpleSeed, DefaultNoScore, void> set5(100,0);
	append(set5,set3);
	SEQAN_TASSERT(length(set5)==2);	
	
	SeedSet<int,SimpleSeed, DefaultNoScore, void> set6(100,0);
	addSeed(set6,0,0,7,Single());
	addSeed(set6,15,15,18,18,Single());
	SEQAN_TASSERT(length(set6)==2);	

	Seed<int,ChainedSeed> seed3(0,0,7);
	SeedSet<int,ChainedSeed, DefaultNoScore, void> set7(100,0);
	addSeed(set7,seed3,Single());
	SEQAN_TASSERT(length(set7)==1);	
	
	SeedSet<int,ChainedSeed, DefaultNoScore, void> set8(100,0);
	addSeedsIt(set8,begin(set7),end(set7),Single());
	addSeed(set8,15,15,6,Single());
	SEQAN_TASSERT(length(set8)==2);

	SeedSet<int,ChainedSeed, DefaultNoScore, void> set9(100,0);
	addSeeds(set9, begin(set8), end(set8), 5, SimpleChain());
	SEQAN_TASSERT(length(set9)==1);

	SEQAN_TASSERT(startDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(endDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(leftDim0(*begin(set9))==0);
	SEQAN_TASSERT(rightDim0(*begin(set9))==20);
	SEQAN_TASSERT(leftDim1(*begin(set9))==0);
	SEQAN_TASSERT(rightDim1(*begin(set9))==20);
	SEQAN_TASSERT(leftDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(rightDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(length(*begin(set9)) == 21);
	

	SeedSet<int,ChainedSeed, DefaultNoScore, void> set10(100,0);
	addSeed(set10,27,28,6,Single());
	Score<int,Simple> matrix(2,-1,-1);
	String<Dna> query =		"AAAAAAAAAAAAAAAAAAAAAACATCGCTTACGCTAT";
	String<Dna> database =	"AAAAAAAAAAAAAAAAAAAAAACAACCTTTAGGCTTT";
	addSeeds(set9, begin(set10), end(set10), matrix, query, database, 4, Chaos());
	
	SEQAN_TASSERT(length(set9)==1);
	SEQAN_TASSERT(startDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(endDiagonal(*begin(set9))==1);
	SEQAN_TASSERT(leftDim0(*begin(set9))==0);
	SEQAN_TASSERT(rightDim0(*begin(set9))==32);
	SEQAN_TASSERT(leftDim1(*begin(set9))==0);
	SEQAN_TASSERT(rightDim1(*begin(set9))==33);
	SEQAN_TASSERT(leftDiagonal(*begin(set9))==1);
	SEQAN_TASSERT(rightDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(length(*begin(set9)) == 33);

	SeedSet<int,ChainedSeed, DefaultNoScore, void> set11(100,0);
	addSeed(set11,0,0,3,Single());

	SeedSet<int,ChainedSeed, DefaultNoScore, void> set12(100,0);
	addSeed(set12,15,15,3,Single());
				//0123456789012345678
	query	 =	 "AAATTTGTTTTGTTTAAA";
	database =	 "AAATTTCTTTTCTTTAAA";
	addSeeds(set11, begin(set12), end(set12), query, database, 12, Blat());
	SEQAN_TASSERT(length(set11)==1);
	Align<String<Dna>, ArrayGaps> aligned;
	getAlignment(*begin(set11),aligned, query, database, matrix);


	clear(set11);
	SEQAN_TASSERT(length(set11)==0);
	addSeed(set11,0,0,3,Single());

	addSeeds(set11, set12, 0, SimpleChain());

	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	
	clear(set11);
	addSeed(set11,0,0,3,Single());
	
	addSeeds(set11, set12, Single());	
	SEQAN_TASSERT(length(set11)==2);

	clear(set11);
	addSeed(set11,0,0,3,Single());

	addSeeds(set11, set12, matrix, query, database, 6, Chaos());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,0,0,3,Single());

	addSeeds(set11, set12, query, database, 9, Blat());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,0,0,7,Single());

	clear(set12);
	addSeed(set12,4,5,7,Single());
	addSeed(set11, 4, 5, 7, 4, Merge());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,0,0,7,Single());
	addSeed(set11, Seed<int,ChainedSeed>(4,5,7), 4, Merge());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	SeedSet<int,SimpleSeed, DefaultNoScore, void> setSimple1(100,0);
	addSeed(setSimple1,0,0,7,Single());
	rightDim0(*begin(setSimple1));
	addSeed(setSimple1, 4, 5, 12, 13, 6, Merge());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==12);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);

	clear(setSimple1);
	addSeed(setSimple1,0,0,7,Single());
	addSeed(setSimple1, 15, 15, 3, 2, SimpleChain());
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==17);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);


	clear(setSimple1);
	addSeed(setSimple1,0,0,7,Single());
	addSeed(setSimple1, 14, 16, 18, 20, 2, SimpleChain());
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==18);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);

	clear(setSimple1);
	addSeed(setSimple1,0,0,7,Single());
	addSeed(setSimple1, Seed<int,SimpleSeed>(15,15,3), 4, SimpleChain());
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==17);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);

	clear(set11);
	addSeed(set11,0,0,7,Single());
	addSeed(set11, 15, 15, 3, matrix, query, database, 5, Chaos());
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,0,0,7,Single());
	addSeed(set11, 15, 15, 3, query,database, 5, Blat());
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,5,5,3,Single());
	extendSeeds(set11, query, database, 2, MatchExtend());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);


	clear(set11);
	addSeed(set11,5,5,3,Single());
	extendSeeds(begin(set11), end(set11), query, database, 2, MatchExtend());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,5,5,3,Single());
	extendSeeds(set11, 1, matrix, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	
	clear(set11);
	addSeed(set11,5,5,3,Single());
	extendSeeds(begin(set11), end(set11), 1, matrix, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	SeedSet<int,ChainedSeed, DefaultNoScore, void> set15(100,0);
	addSeed(set15,5,5,3,Single());
	extendSeeds(begin(set15), end(set15), 1, matrix, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(rightDim0(*begin(set15))==10);
	SEQAN_TASSERT(leftDim0(*begin(set15))==0);
}


void Test_SeedSet_seedScore()
{
	Score<int,Simple> matrix(2,-1,-1, 0);
	SeedSet<int,SimpleSeed,DefaultScore, void> set1;
	setScoreMatrix(set1,matrix);
	matrix = getScoreMatrix(set1);
	SeedSet<int,ChainedSeed,DefaultScore, void> set2;
	SeedSet<int,SimpleSeed,DefaultScore, void> set3(100,0,matrix);
	SeedSet<int,ChainedSeed,DefaultScore, void> set4(100,0,matrix);
	
	Iterator<SeedSet<int,ChainedSeed,DefaultScore, void>, Standard>::Type it_multi;
	Iterator<SeedSet<int,SimpleSeed,DefaultScore, void>, Standard>::Type it_simple;
	SEQAN_TASSERT(qualityValue(set3)==0);
	SEQAN_TASSERT(maximumDistance(set3)==100);
	setQualityValue(set3,1);
	setMaximumDistance(set3,99);
	SEQAN_TASSERT(qualityValue(set3)==1);
	SEQAN_TASSERT(maximumDistance(set3)==99);
	SEQAN_TASSERT(length(set3)==0);	
	//SEQAN_TASSERT(capacity(set3)==1);

	Seed<int,SimpleSeed> seed(0,0,3);
	Seed<int,SimpleSeed> seed2(3,4,3);
	appendValue(set3,seed,3);
	appendValue(set3,seed2,4);
    it_simple = begin(set3);
	SEQAN_TASSERT(seedScore(it_simple)==3);	
	SEQAN_TASSERT(seedScore(++it_simple)==4);	

	String<int> seedScores;
	appendValue(seedScores,4);
	appendValue(seedScores,5);
	appendValue(seedScores,3);
	SeedSet<int,SimpleSeed,DefaultScore, void> set5(100,0,matrix);

	append(set5,set3, seedScores);
	SEQAN_TASSERT(length(set5)==2);	
        it_simple = begin(set5);
	SEQAN_TASSERT(seedScore(it_simple)==4);	
	SEQAN_TASSERT(seedScore(++it_simple)==5);	

	SeedSet<int,SimpleSeed,DefaultScore, void> set6(100,0,matrix);
	addSeed(set6,0,0,7,Single());
	addSeed(set6,15,15,18,18,6,Single());
	SEQAN_TASSERT(length(set6)==2);	
        it_simple = begin(set6);
	SEQAN_TASSERT(seedScore(it_simple)==14);	
	SEQAN_TASSERT(seedScore(++it_simple)==6);	

	clear(set6);

	Seed<int,ChainedSeed> seed3(0,0,7);
	SeedSet<int,ChainedSeed,DefaultScore, void> set7(100,0,matrix);
	addSeed(set7,seed3,7,Single());
	SEQAN_TASSERT(length(set7)==1);	
	SEQAN_TASSERT(seedScore(begin(set7))==7);
	
	SeedSet<int,ChainedSeed,DefaultScore, void> set8(100,0,matrix);
	addSeeds(set8,begin(set7),end(set7), begin(seedScores), Single());
	addSeed(set8,10,10,6,7,Single());
	SEQAN_TASSERT(length(set8)==2);
	it_multi = begin(set8);
	SEQAN_TASSERT(seedScore(it_multi)== 4);
	SEQAN_TASSERT(rightDim0(*it_multi)==6);
    SEQAN_TASSERT(seedScore(++it_multi)==7);
	SEQAN_TASSERT(rightDim0(*it_multi)==15);

	clear(set8);
	addSeedSet(set8,set7, Single());
	addSeed(set8,10,10,6,7,Single());
	SEQAN_TASSERT(length(set8)==2);
        it_multi = begin(set8);
	SEQAN_TASSERT(seedScore(it_multi)== 7);
	SEQAN_TASSERT(rightDim0(*it_multi)==6);
        SEQAN_TASSERT(seedScore(++it_multi)==7);
	SEQAN_TASSERT(rightDim0(*it_multi)==15);

	SeedSet<int,ChainedSeed,DefaultScore, void> set9(100,0,matrix);
	addSeeds(set9, begin(set8), end(set8), begin(seedScores), 5, SimpleChain());
	SEQAN_TASSERT(length(set9)==1);
	SEQAN_TASSERT(seedScore(begin(set9))==3);

	SEQAN_TASSERT(startDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(endDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(leftDim0(*begin(set9))==0);
	SEQAN_TASSERT(rightDim0(*begin(set9))==15);
	SEQAN_TASSERT(leftDim1(*begin(set9))==0);
	SEQAN_TASSERT(rightDim1(*begin(set9))==15);
	SEQAN_TASSERT(leftDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(rightDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(length(*begin(set9)) == 16);
	
	SeedSet<int,ChainedSeed,DefaultScore, void> set10(100,0,matrix);
	addSeed(set10,27,28,6,5,Single());
	
			       //0123456789012345678901234567890123456
	String<Dna> query =	"AAAAAAAAAAAAAAAAAAAAAACATCGCTTACGCTAT";
	String<Dna> database =	"AAAAAAAAAAAAAAAAAAAAAACAACCTTAGGCTTT";

	addSeeds(set9, begin(set10), end(set10), begin(seedScores), query, database, 4, Chaos());
	SEQAN_TASSERT(seedScore(begin(set9))==22);
	SEQAN_TASSERT(length(set9)==1);
	SEQAN_TASSERT(startDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(endDiagonal(*begin(set9))==1);
	SEQAN_TASSERT(leftDim0(*begin(set9))==0);
	SEQAN_TASSERT(rightDim0(*begin(set9))==32);
	SEQAN_TASSERT(leftDim1(*begin(set9))==0);
	SEQAN_TASSERT(rightDim1(*begin(set9))==33);
	SEQAN_TASSERT(leftDiagonal(*begin(set9))==1);
	SEQAN_TASSERT(rightDiagonal(*begin(set9))==0);
	SEQAN_TASSERT(length(*begin(set9)) == 33);

	SeedSet<int,ChainedSeed,DefaultScore, void> set11(100,0,matrix);
	addSeed(set11,0,0,3,5,Single());

	SeedSet<int,ChainedSeed,DefaultScore, void> set12(100,0,matrix);
	addSeed(set12,15,15,3,5,Single());

	SeedSet<int,SimpleSeed,DefaultScore, void> set11a(100,0,matrix);
	addSeed(set11a,0,0,3,5,Single());

	SeedSet<int,SimpleSeed,DefaultScore, void> set12a(100,0,matrix);
	addSeed(set12a,15,15,3,5,Single());

	query	 =	 "AAATTTGTTTTGTTTAAA";
	database =	 "AAATTTCTTTTCTTTAAA";
	addSeeds(set11, begin(set12), end(set12), begin(seedScores), query, database, 12, Blat());
	SEQAN_TASSERT(length(set11)==1);
	Align<String<Dna>, ArrayGaps> aligned;

	SEQAN_TASSERT(seedScore(begin(set11))==25);

	String<Seed<int,ChainedSeed> > seedString;
	appendValue(seedString,Seed<int, ChainedSeed>(15,15,3));

	clear(set11);
	SEQAN_TASSERT(length(set11)==0);
	addSeed(set11,0,0,3,5,Single());
	addSeeds(set11, seedString, seedScores, 3, SimpleChain());
	SEQAN_TASSERT(seedScore(begin(set11))==-15);
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	SEQAN_TASSERT(length(set11)==0);
	addSeed(set11,0,0,3,Single());
	addSeed(set11, 15, 15, 3, 3, SimpleChain());
	SEQAN_TASSERT(seedScore(begin(set11))==-12);
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,0,0,3,5,Single());
	addSeedSet(set11, set12, 3, SimpleChain());
	SEQAN_TASSERT(seedScore(begin(set11))==-14);
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,0,0,3,6,Single());
	addSeeds(set11, set12, seedScores,  Single());	
	SEQAN_TASSERT(length(set11)==2);

	clear(set11);
	addSeed(set11,0,0,3,5,Single());
	addSeeds(set11, set12, seedScores, query, database, 6, Chaos());	
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==27);
	
	clear(set11a);
	addSeed(set11a,0,0,3,5,Single());
	addSeeds(set11a, set12a, seedScores, query, database, 6, Chaos());	
	SEQAN_TASSERT(length(set11a)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11a))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11a))==0);
	SEQAN_TASSERT(seedScore(begin(set11a))==27);

	clear(set11a);
	addSeed(set11a,0,0,3,5,Single());
	addSeed(set11a,15,15,3,4,query, database, 5, Chaos());	
	SEQAN_TASSERT(length(set11a)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11a))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11a))==0);
	SEQAN_TASSERT(seedScore(begin(set11a))==27);
	
	clear(set11);
	addSeed(set11,0,0,3,Single());
	addSeed(set11, 15, 15, 3, query, database, 6, Chaos());	
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==30);

	clear(set11a);
	addSeed(set11a,0,0,3,Single());
	addSeed(set11a, 15, 15, 3, query, database, 6, Chaos());	
	SEQAN_TASSERT(length(set11a)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11a))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11a))==0);
	SEQAN_TASSERT(seedScore(begin(set11a))==30);

	clear(set11);
	addSeed(set11,0,0,3,5,Single());
	addSeedSet(set11, set12, query, database, 6, Chaos());	
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==28);

	clear(set11);
	addSeed(set11,0,0,3,5,Single());
	addSeeds(set11, set12, seedScores, query, database, 9, Blat());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==25);

	clear(set11);
	addSeed(set11,0,0,3,Single());
	addSeed(set11, 15, 15, 3, query, database, 9, Blat());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==28);

	clear(set11a);
	addSeed(set11a,0,0,3,Single());
	addSeed(set11a, 15, 15, 3, query, database, 9, Blat());
	SEQAN_TASSERT(length(set11a)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11a))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11a))==0);
	SEQAN_TASSERT(seedScore(begin(set11a))==28);

	clear(set11);
	addSeed(set11,0,0,3,5,Single());
	addSeedSet(set11, set12, query, database, 9, Blat());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==26);

	clear(set11a);
	addSeed(set11a,0,0,3,5,Single());
	addSeedSet(set11a, set12a, query, database, 9, Blat());
	SEQAN_TASSERT(length(set11a)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11a))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11a))==0);
	SEQAN_TASSERT(seedScore(begin(set11a))==26);


	clear(set11);
	addSeed(set11,0,0,7,14,Single());

	clear(set12);
	addSeed(set12,4,5,7,14,Single());
	
	addSeed(set11, 4, 5, 7, 14, 4, Merge());
	SEQAN_TASSERT(length(set11)==1);
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==21);

	clear(set11);
	addSeed(set11,0,0,7,14,Single());
	addSeed(set11,4,5,7,14,2,Merge());
	addSeed(set11, Seed<int,ChainedSeed>(3,4,10), 20, 5,Merge());
	SEQAN_TASSERT(length(set11)==1);
        it_multi = begin(set11);
	//SEQAN_TASSERT(rightDim0(*++it_multi)==12);
	SEQAN_TASSERT(leftDim0(*it_multi)==0);
	SEQAN_TASSERT(seedScore(it_multi)==25);


	clear(set11);
	addSeed(set11,0,0,7,Single());
	addSeed(set11,4,5,7,2,Merge());
	addSeed(set11, Seed<int,ChainedSeed>(3,4,10), 20, 5,Merge());
	SEQAN_TASSERT(length(set11)==1);
        it_multi =begin(set11);
	SEQAN_TASSERT(leftDim0(*it_multi)==0);
	SEQAN_TASSERT(seedScore(it_multi)==25);


	clear(set11);
	addSeed(set11,0,0,7,5,Single());
	addSeed(set11, 15, 15, 3, 2,5, SimpleChain());
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	

	SeedSet<int,SimpleSeed,DefaultScore, void> setSimple1(100,0,matrix);
	addSeed(setSimple1,0,0,7,14,Single());
	rightDim0(*begin(setSimple1));
	addSeed(setSimple1, 5, 4, 12, 13, 14, 5,Merge());
	SEQAN_TASSERT(length(setSimple1)==1);
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==12);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);
	SEQAN_TASSERT(seedScore(begin(setSimple1))==21);

	clear(setSimple1);
	addSeed(setSimple1,0,0,7,5,Single());
	addSeed(setSimple1, 15, 15, 3, 2,5, SimpleChain());
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==17);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);

	clear(setSimple1);
	addSeed(setSimple1,0,0,7,5, Single());
	addSeed(setSimple1, 14, 16, 18, 20, 2, 5, SimpleChain());
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==18);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);

	clear(setSimple1);
	addSeed(setSimple1,0,0,7,5, Single());
	addSeed(setSimple1, Seed<int,SimpleSeed>(15,15,3), 4, 5, SimpleChain());
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==17);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);


	clear(set11);
	addSeed(set11,0,0,7,5, Single());
	addSeed(set11, 15, 15, 3, 5, query, database, 5, Chaos());
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,0,0,7,5, Single());
	addSeed(set11, 15, 15, 3, 5, query, database, 5, Blat());
	SEQAN_TASSERT(rightDim0(*begin(set11))==17);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	clear(set11);
	addSeed(set11,5,5,3,5, Single());
	extendSeedsScore(set11, seedScores, matrix, query, database, 2, MatchExtend());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScores[0]==20);
	seedScores[0] = 4;

	clear(set11);
	addSeed(set11,5,5,3,5, Single());
	extendSeedsScore(set11, query, database, 2, MatchExtend());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==21);

	clear(set11);
	addSeed(set11,5,5,3,5, Single());
	extendSeedsScore(begin(set11), end(set11), begin(seedScores), matrix, query, database, 2, MatchExtend());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	
	seedScores[0] = 4;
	clear(set11);
	addSeed(set11,5,5,3,5,Single());
	extendSeedsScore(set11, seedScores, 1, matrix, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScores[0]==20);
	
	clear(set11);
	addSeed(set11,5,5,3,5,Single());
	extendSeedsScore(set11,1, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);
	SEQAN_TASSERT(seedScore(begin(set11))==21);

	clear(setSimple1);
	addSeed(setSimple1,5,5,3,5,Single());
	extendSeedsScore(setSimple1,1, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(rightDim0(*begin(setSimple1))==10);
	SEQAN_TASSERT(leftDim0(*begin(setSimple1))==0);
	SEQAN_TASSERT(seedScore(begin(setSimple1))==21);

	clear(set11);
	addSeed(set11,5,5,3,5,Single());
	extendSeedsScore(begin(set11), end(set11), begin(seedScores), 1, matrix, query, database, 2, UngappedXDrop());
	SEQAN_TASSERT(rightDim0(*begin(set11))==10);
	SEQAN_TASSERT(leftDim0(*begin(set11))==0);

	seedScores[0] = 5;
	SeedSet<int,SimpleSeed,DefaultScore, void> setSimple2(100,0,matrix);
	addSeed(setSimple2,5,5,3,5,Single());
	extendSeedsScore(begin(setSimple2), end(setSimple2), begin(seedScores), 1, matrix, query, database, 2, GappedXDrop());
	SEQAN_TASSERT(seedScores[0]==20);

	seedScores[0] = 5;
	SeedSet<int,ChainedSeed,DefaultScore, void> setSimple3(100,0,matrix);
	addSeed(setSimple3,5,5,3,5,Single());
	extendSeedsScore(begin(setSimple3), end(setSimple3), begin(seedScores), 1, matrix, query, database, 2, GappedXDrop());
	SEQAN_TASSERT(seedScores[0]==20);
}

void Main_SeedSet(){
	SEQAN_TREPORT("TEST BEGIN")
	Test_SeedSet_base();
	Test_SeedSet_seedScore();
	debug::verifyCheckpoints("projects/library/seqan/seeds/seedSet_base.h");
	debug::verifyCheckpoints("projects/library/seqan/seeds/seedSet_score.h");
	SEQAN_TREPORT("TEST END")
}
