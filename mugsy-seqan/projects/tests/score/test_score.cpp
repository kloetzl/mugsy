#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/score.h>

using namespace std;
using namespace seqan;

#define TEST_PATH "projects/tests/score/"
#define LIB_PATH "projects/library/seqan/score/"

//////////////////////////////////////////////////////////////////////////////

template <typename TScore1, typename TScore2>
void testCompareAAMatrices(TScore1 const & mat1, 
						   TScore2 const & mat2)
{
	AminoAcid a, b;
	for (a = 'A'; a <= '*'; ++a)
	{
		for (b = 'A'; b <= '*'; ++b)
		{
			SEQAN_TASSERT(score(mat1, a, b) == score(mat2, a, b))
		}
	}
}

void testMatrixScore()
{
	Score<int, ScoreMatrix<> > sc;
	String<char> meta;
	loadScoreMatrix(sc, TEST_PATH "BLOSUM62", meta);

	//compare the matrix loaded from file with the build-in matrix
	Blosum62 blosum62;
	testCompareAAMatrices(blosum62, sc);

	SEQAN_TASSERT(scoreGapExtend(blosum62) == -1);
	SEQAN_TASSERT(scoreGapOpen(blosum62) == scoreGapExtend(blosum62));
	SEQAN_TASSERT(scoreGap(blosum62) == scoreGapExtend(blosum62));

	//store and load again
	FILE * fl = fopen(TEST_PATH "testfile.txt", "wb");
	write(fl, sc, meta);
	fclose(fl);

	Score<int, ScoreMatrix<> > sc2;
	String<char> meta2;
	loadScoreMatrix(sc2, TEST_PATH "testfile.txt", meta2);
	testCompareAAMatrices(sc, sc2);
	SEQAN_TASSERT(meta == meta2)

	//store and load again build-in matrix
	fl = fopen(TEST_PATH "testfile.txt", "wb");
	write(fl, Blosum62());
	fclose(fl);

	loadScoreMatrix(sc2, TEST_PATH "testfile.txt");
	testCompareAAMatrices(sc, Blosum62());

	//test setScore()
	setScore(sc, 'A', '*', 100);
	SEQAN_TASSERT(score(sc, 'A', '*') == 100)
}

//////////////////////////////////////////////////////////////////////////////
//compare to http://www.bioinformatics.nl/tools/pam.html

void testScorePAM()
{
//	PamJones pam;

	PamDayhoff pam;
	//Score<int, Pam<> > pam(250, -1, 0);
	//Score<int, Pam<AminoAcid, Pam_Data_Dayhoff_MDM78> > pam;

	SEQAN_TASSERT(getDist(pam) == 250);
	SEQAN_TASSERT(scoreGapExtend(pam) == -1);
	SEQAN_TASSERT(scoreGapOpen(pam) == scoreGapExtend(pam));
	SEQAN_TASSERT(scoreGap(pam) == scoreGapExtend(pam));

	//store and load again build-in matrix
	FILE * fl = fopen(TEST_PATH "testfile.txt", "wb");
	write(fl, pam);
	fclose(fl);

	Score<int, ScoreMatrix<> > sc;
	loadScoreMatrix(sc, TEST_PATH "testfile.txt");
	testCompareAAMatrices(sc, PamDayhoff());
}


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	testScorePAM();

	testMatrixScore();

//____________________________________________________________________________

	SEQAN_TREPORT("TEST END")

	return 0;
}
