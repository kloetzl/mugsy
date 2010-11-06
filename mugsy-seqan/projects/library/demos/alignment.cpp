///A tutorial about global alignments.
#include <iostream>
#include <seqan/align.h>
#include <seqan/graph_align.h>

using namespace seqan;

int main()
{
///Two DNA sequences that shall be aligned.
	typedef String<Dna> TSequence;
	TSequence seq1 = "atcgaatgcgga";
	TSequence seq2 = "actcgttgca";
///Scoring objects are used to define a scoring scheme.
///In this case, affine gap costs with match = 0, mismatch = -1, gapextend = -1 and gapopen = -2.
	Score<int> score(0, -1, -1, -2);
///Example 1: We use @Class.Align@ to align the two sequences. 
///Since we do not specify an @Tag.Global Alignment Algorithms|algorithm tag@ when we call @Function.globalAlignment@, 
///a suitable algorithm (@Tag.Global Alignment Algorithms|Gotoh@) is automatically choosen.
	Align<TSequence, ArrayGaps> align;
	resize(rows(align), 2);
	assignSource(row(align, 0), seq1);
	assignSource(row(align, 1), seq2);

	::std::cout << "Score = " << globalAlignment(align, score) << ::std::endl;
	::std::cout << align << ::std::endl;
///Example 2: We now choose explicitely the algorithm @Tag.Global Alignment Algorithms|MyersHirschberg@.
///Since this algorithm always works on Levenshtein distance, $score$ is ignored here.
///Therefore, this algorithm computes a different alignment and returns a different score.
	::std::cout << "Score = " << globalAlignment(align, score, MyersHirschberg()) << ::std::endl;
	::std::cout << align << ::std::endl;
///Example 3: We now do the same as in case 1, but now we use an @Spec.Alignment Graph@ for storing the alignment.
///Here we use @Tag.Global Alignment Algorithms|Gotoh's algorithm@.
	typedef StringSet<TSequence, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

	TStringSet string_set;
	appendValue(string_set, seq1);
	appendValue(string_set, seq2);
	TAlignmentGraph alignment_graph(string_set);

	::std::cout << "Score = " << globalAlignment(alignment_graph, score, Gotoh()) << ::std::endl;
	::std::cout << alignment_graph << ::std::endl;
	return 0;
}
