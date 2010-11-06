#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>

using namespace std;
using namespace seqan;

template <typename TAlphabet>
String<TAlphabet> generate_random(int length_of_sequence)
{
	// init string
	String<TAlphabet> ret;
	resize(ret,length_of_sequence);
	int alphabet_size = ValueSize<TAlphabet>::VALUE;
	// generate random sequence of length "length_of_sequence"
	for (int i = 0; i < length_of_sequence; ++i)
		ret[i] =  static_cast<int>((alphabet_size * static_cast<int>(rand()) / (static_cast<int>(RAND_MAX) + 1)));

	return ret;
}

template <typename TAlphabet>
String<TAlphabet> generate_second_sequence(int error_count,String<TAlphabet> copy_of_org)
{
	int length_of_org = length(copy_of_org);
	
	int alphabet_size = ValueSize<TAlphabet>::VALUE;

	for(int i = 0;i < error_count;++i)
	{
		// introduce errors into sequence
		// 1. choose position 
		int pos = static_cast<int>((int)length_of_org * rand() / (RAND_MAX +  1.0));
		
		// generate new char
		TAlphabet new_char = static_cast<TAlphabet>((alphabet_size * static_cast<int>(rand()) / (static_cast<int>(RAND_MAX) + 1)));
		
		// replace char
		copy_of_org[pos] = new_char;
	}
	return copy_of_org;
}


template <typename TAlphabet>
void erase_sequence_parts(int erase_count,String<TAlphabet> & sequence)
{
	// erase single characters
	int len = length(sequence) - 1;
	for(int i = 0;i < erase_count;++i)
	{
		// calc position
		int pos = static_cast<int>((int)(len - i) * rand() / (RAND_MAX +  1.0));
		erase(sequence,pos,pos+1);
	}
}

#define ALPHABET Dna

int edit_distance(Align<String<ALPHABET>, ArrayGaps> & ali)
{
	int len_ali = length(row(ali,0));

	Iterator<Row<Align<String<ALPHABET>, ArrayGaps> >::Type >::Type ali_row_0 = iter(row(ali, 0), 0);
	Iterator<Row<Align<String<ALPHABET>, ArrayGaps> >::Type >::Type ali_row_1 = iter(row(ali, 1), 0);

	int score = 0;
	int i;
	// iteration ueber das alignment
	for(i = 0;i < len_ali;++i)
	{
		if(isGap(ali_row_0))
		{
			--score;
		}
		else if(isGap(ali_row_1))
		{
			--score;
		}
		else if(value(ali_row_0) != value(ali_row_1))
		{
			--score;
		}

		goNext(ali_row_0);
		goNext(ali_row_1);
	}

	return score;
}


// [Test	]
void Main_TestMyers()
{
	SEQAN_TREPORT("TEST MYERS ALIGN BEGIN")

	int nw_score,m_score,hm_score;
	int test_repeat = 15;
	int test_count = 0;

	// ---------------------- short sequences ----------------------
	SEQAN_TREPORT("TEST MYERS SCORE AND ALIGN WITH SHORT SEQUENCES -- BEGIN")
	printf("\n");
	while(test_count < test_repeat)
	{
		printf("\rrun %i of %i",(test_count+1),test_repeat);

		// create random sequences 
		String<ALPHABET> s_str0 = generate_random<ALPHABET>(20);
		String<ALPHABET> s_str1 = generate_second_sequence<ALPHABET>(3,s_str0);
		erase_sequence_parts(3,s_str1);

		// test alignment with random sequences 
		// use needleman wunsch as reference
		Align<String<ALPHABET>, ArrayGaps> s_nw_ali;
		resize(rows(s_nw_ali), 2);
		assignSource(row(s_nw_ali, 0), s_str0);
		assignSource(row(s_nw_ali, 1), s_str1);

		nw_score = globalAlignment(s_nw_ali,SimpleScore());

		// compute only score
		Align<String<ALPHABET>, ArrayGaps> s_m_ali;
		resize(rows(s_m_ali), 2);
		assignSource(row(s_m_ali, 0), s_str0);
		assignSource(row(s_m_ali, 1), s_str1);

		m_score = globalAlignment(s_m_ali,SimpleScore(), MyersBitVector());
		SEQAN_TASSERT(nw_score == m_score);
		
		// compute complete alignments
		Align<String<ALPHABET>, ArrayGaps> s_hm_ali;
		resize(rows(s_hm_ali), 2);
		assignSource(row(s_hm_ali, 0), s_str0);
		assignSource(row(s_hm_ali, 1), s_str1);

		hm_score = globalAlignment(s_hm_ali,SimpleScore(), MyersHirschberg());
		SEQAN_TASSERT(nw_score == hm_score);
		SEQAN_TASSERT(edit_distance(s_hm_ali) == hm_score);

		++test_count;
	}
	printf("\n\n");
	SEQAN_TREPORT("TEST MYERS SCORE AND ALIGN WITH SHORT SEQUENCES -- END")

	// ---------------------- long sequences ----------------------

	SEQAN_TREPORT("TEST MYERS SCORE AND ALIGN WITH LONG SEQUENCES -- BEGIN")
	printf("\n");
	test_count = 0;

	while(test_count < test_repeat)
	{
	
		printf("\rrun %i of %i",(test_count+1),test_repeat);

		// create random sequences 
		String<ALPHABET> l_str0 = generate_random<ALPHABET>(200);
		String<ALPHABET> l_str1 = generate_second_sequence<ALPHABET>(10,l_str0);
		erase_sequence_parts(5,l_str1);
		
		// test alignment with random sequences 
		// use needleman wunsch as reference
		Align<String<ALPHABET>, ArrayGaps> l_nw_ali;
		resize(rows(l_nw_ali), 2);
		assignSource(row(l_nw_ali, 0), l_str0);
		assignSource(row(l_nw_ali, 1), l_str1);

		nw_score = globalAlignment(l_nw_ali,SimpleScore());

		// compute only score
		Align<String<ALPHABET>, ArrayGaps> l_m_ali;
		resize(rows(l_m_ali), 2);
		assignSource(row(l_m_ali, 0), l_str0);
		assignSource(row(l_m_ali, 1), l_str1);

		m_score = globalAlignment(l_m_ali,SimpleScore(), MyersBitVector());

		SEQAN_TASSERT(nw_score == m_score);
		
		// compute complete alignments
		Align<String<ALPHABET>, ArrayGaps> l_hm_ali;
		resize(rows(l_hm_ali), 2);
		assignSource(row(l_hm_ali, 0), l_str0);
		assignSource(row(l_hm_ali, 1), l_str1);

		hm_score = globalAlignment(l_hm_ali, SimpleScore(), MyersHirschberg());

		SEQAN_TASSERT(nw_score == hm_score);
		SEQAN_TASSERT(edit_distance(l_hm_ali) == hm_score);

		++test_count;
	}
	printf("\n\n");
	SEQAN_TREPORT("TEST MYERS SCORE AND ALIGN WITH LONG SEQUENCES -- END")

	SEQAN_TREPORT("TEST MYERS ALIGN END")

	// ---------------------- hirschberg algorithm ----------------------
	SEQAN_TREPORT("TEST HIRSCHBERG -- BEGIN")
	test_count = 0;

	while(test_count < test_repeat)
	{
	
		printf("\rrun %i of %i",(test_count+1),test_repeat);

		// create random sequences 
		String<ALPHABET> str0 = generate_random<ALPHABET>(20);
		String<ALPHABET> str1 = generate_second_sequence<ALPHABET>(2,str0);

		erase_sequence_parts(5,str1);
		
		// test alignment with random sequences 
		// use needleman wunsch as reference
		Align<String<ALPHABET>, ArrayGaps> nw_ali;
		resize(rows(nw_ali), 2);
		assignSource(row(nw_ali, 0), str0);
		assignSource(row(nw_ali, 1), str1);

		nw_score = globalAlignment(nw_ali,SimpleScore());
		
		// compute complete alignments with hirschberg algorithm
		Align<String<ALPHABET>, ArrayGaps> hirsch_ali;
		resize(rows(hirsch_ali), 2);
		assignSource(row(hirsch_ali, 0), str0);
		assignSource(row(hirsch_ali, 1), str1);

		hm_score = globalAlignment(hirsch_ali,SimpleScore(), Hirschberg());

		SEQAN_TASSERT(nw_score == hm_score);
		SEQAN_TASSERT(edit_distance(hirsch_ali) == hm_score);

		++test_count;
	}
	printf("\n");
	SEQAN_TREPORT("TEST HIRSCHBERG -- END")

	debug::verifyCheckpoints("projects/library/seqan/align/align_myers.h");
	debug::verifyCheckpoints("projects/library/seqan/align/align_hirschberg.h");
	debug::verifyCheckpoints("projects/library/seqan/align/hirschberg_set.h");
}

