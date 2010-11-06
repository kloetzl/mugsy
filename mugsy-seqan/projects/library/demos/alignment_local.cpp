///A tutorial about local alignments.
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main()
{
///Example 1: This program applies the Smith-Waterman algorithm to compute the best local alignment between two given sequences.
	Align< String<char> > ali;
	appendValue(rows(ali), "aphilologicaltheorem");
	appendValue(rows(ali), "bizarreamphibology");
    ::std::cout << "Score = " << localAlignment(ali, Score<int>(3,-3,-2, -2), SmithWaterman()) << ::std::endl;
	::std::cout << ali;
	::std::cout << "Aligns Seq1[" << sourceBeginPosition(row(ali, 0)) << ":" << (sourceEndPosition(row(ali, 0))-1) << "]";
	::std::cout << " and Seq2[" << sourceBeginPosition(row(ali, 1)) << ":" <<  (sourceEndPosition(row(ali, 1))-1) << "]" << ::std::endl << ::std::endl;


///Example 2: This program applies the Waterman-Eggert algorithm to compute all non-overlapping local alignments with score better than 2.
	Align< String<Dna> > ali2;
	appendValue(rows(ali2), "ataagcgtctcg");
	appendValue(rows(ali2), "tcatagagttgc");

	LocalAlignmentFinder<> finder(ali2);
	Score<int> scoring(2, -1, -2, 0);
	while (localAlignment(ali2, finder, scoring, 2, WatermanEggert())) {
		::std::cout << "Score = " << getScore(finder) << ::std::endl;
		::std::cout << ali2;
		::std::cout << "Aligns Seq1[" << sourceBeginPosition(row(ali2, 0)) << ":" << (sourceEndPosition(row(ali2, 0))-1) << "]";
		::std::cout << " and Seq2[" << sourceBeginPosition(row(ali2, 1)) << ":" <<  (sourceEndPosition(row(ali2, 1))-1) << "]" << ::std::endl << ::std::endl;
	}

	return 0;
}
