#include <iostream>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/find_motif.h>

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template<typename TIter, typename TString, typename TType>
bool
isOOPSMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) & (((int)(ds_end-ds_iter))==t) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<2) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd == (int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<=(int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t==0)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

template<typename TIter, typename TString, typename TType>
bool
isOMOPSMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) & (((int)(ds_end-ds_iter))==t) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<1) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd== (int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<=(int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t==0)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

template<typename TIter, typename TString, typename TType>
bool
isZOOPSMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int upper_limit = t-((int) floor(t*(ZOOPS().threshold)+0.5)-1);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<2) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd== (int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<= (int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t<=upper_limit)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

template<typename TIter, typename TString, typename TType>
bool
isTCMMotif(TIter ds_iter, TIter ds_end, TString const & motif, TType const & d, bool & is_exact)
{
	bool isRelevantMotif = false;
	int t = (int) (ds_end-ds_iter);
	int upper_limit = t-((int) floor(t*(TCM().threshold)+0.5)-1);
	int counter = 0;
	int hd = 0;
	typename Size<TString>::Type l = length(motif);
	typename Iterator<TString>::Type seq_iter, seq_end, l_mer_begin;
	while( (ds_iter!=ds_end) )
	{
		seq_iter = begin(*ds_iter);
		seq_end = seq_iter+length(*ds_iter)-l+1;
		counter = 0;
		while( (seq_iter!=seq_end) & (counter<1) )
		{
			l_mer_begin = begin(motif);
			hd = hammingDistance<int>(seq_iter, seq_iter+l, l_mer_begin);
			if( is_exact & (hd==(int) d) )
			{
				++counter;
			}
			else if( !is_exact & (hd<=(int) d))
			{
				++counter;
			}
			++seq_iter;
		}
		if(counter==1)
		{
			--t;
		}
		++ds_iter;
	}
	if(t<=upper_limit)
	{
		isRelevantMotif = true;
	}
	return isRelevantMotif;
}

//////////////////////////////////////////////////////////////////////////////

void Test_exactAlgorithms()
{
	//Testing PMS1 & PMSP algorithm

	unsigned int t = 0;      //number of sequences
	unsigned int l = 0;		//length of motif
	unsigned int d = 0;		//number of substitutions
	bool is_exact = false;	//size of Hamming distance
	unsigned int i = 0;

//____________________________________________________________________________
// Test1 - Search for OOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	t = 3;
	l = 4;		
	d = 1;		
	is_exact = true;	

	String<DnaString> dataset1;
	resize(dataset1, t);
	dataset1[0] = "ACAGCA";
	dataset1[1] = "AGGCAG";
	dataset1[2] = "TCAGTC";

	//Application of PMS1-OOPS
	MotifFinder<Dna, PMS1> motif_finder1(l,d,is_exact);
	findMotif(motif_finder1,dataset1,OOPS());

	//Application of PMSP-OOPS
	MotifFinder<Dna, PMSP> motif_finder2(l,d,is_exact);
	findMotif(motif_finder2,dataset1,OOPS());

	SEQAN_TASSERT(length(motif_finder1.set_of_motifs)==length(motif_finder2.set_of_motifs));
	for(i=0; i<length(motif_finder1.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(motif_finder1.set_of_motifs[i]==motif_finder2.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test2 - Search for OMOPS motifs on a small set of nucleotide sequences
//         given the inexact Hamming distance (<=d)

	l = 6;		
	d = 2;		
	is_exact = false;	

	String<DnaString> dataset2;
	appendValue(dataset2,DnaString("GCTGGACGTG"));
	appendValue(dataset2,DnaString("TCTAGACATA"));
	appendValue(dataset2,DnaString("AGTGGGGGAC"));
	appendValue(dataset2,DnaString("CTAGTCAAGA"));
	appendValue(dataset2,DnaString("CTCGAGGGGT"));

	//Application of PMS1-OMOPS
	MotifFinder<Dna, PMS1> motif_finder3(l,d,is_exact);
	findMotif(motif_finder3,dataset2,OMOPS());

	//Application of PMSP-OMOPS
	MotifFinder<Dna, PMSP> motif_finder4(l,d,is_exact);
	findMotif(motif_finder4,dataset2,OMOPS());

	SEQAN_TASSERT(length(motif_finder3.set_of_motifs)==length(motif_finder4.set_of_motifs));
	for(i=0; i<length(motif_finder3.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(motif_finder3.set_of_motifs[i]==motif_finder4.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test3 - Search for ZOOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	l = 6;		
	d = 1;		
	is_exact = true;	

	String<DnaString> dataset3;
	appendValue(dataset3,DnaString("AGCCGTCTGA"));
	appendValue(dataset3,DnaString("TCCAGGCAAG"));
	appendValue(dataset3,DnaString("GAACGTCCAA"));
	appendValue(dataset3,DnaString("GCTTTCTAAC"));
	appendValue(dataset3,DnaString("AGTAGCTCGC"));

	//Application of PMS1-ZOOPS
	MotifFinder<Dna, PMS1> motif_finder5(l,d,is_exact);
	findMotif(motif_finder5,dataset3,ZOOPS());

	//Application of PMSP-ZOOPS
	MotifFinder<Dna, PMSP> motif_finder6(l,d,is_exact);
	findMotif(motif_finder6,dataset3,ZOOPS());

	SEQAN_TASSERT(length(motif_finder5.set_of_motifs)==length(motif_finder6.set_of_motifs));
	for(i=0; i<length(motif_finder5.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(motif_finder5.set_of_motifs[i]==motif_finder6.set_of_motifs[i]);
	}

//____________________________________________________________________________
// Test4 - Search for TCM motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	l = 5;		
	d = 1;		

	String<DnaString> dataset4;
	appendValue(dataset4,DnaString("CAATTAACTC"));
	appendValue(dataset4,DnaString("ATAAACAGTG"));
	appendValue(dataset4,DnaString("GAATGCATTG"));

	//Application of PMS1-TCM
	MotifFinder<Dna, PMS1> motif_finder7(l,d,is_exact);
	findMotif(motif_finder7,dataset4,TCM());

	//Application of PMSP-TCM
	MotifFinder<Dna, PMSP> motif_finder8(l,d,is_exact);
	findMotif(motif_finder8,dataset4,TCM());
	
	SEQAN_TASSERT(length(motif_finder7.set_of_motifs)==length(motif_finder8.set_of_motifs));
	for(i=0; i<length(motif_finder7.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(motif_finder7.set_of_motifs[i]==motif_finder8.set_of_motifs[i]);
	}

}

void Test_approximationAlgorithms()
{
	//Testing Projection & ePatternBranching algorithm
	
	unsigned int t = 0;      //number of sequences
	unsigned int n = 0;		//length of sequence
	unsigned int l = 0;		//length of motif
	unsigned int d = 0;		//number of substitutions
	bool is_exact = false;	//size of Hamming distance
	unsigned int m =0;		//total number of possible l-mers
	unsigned int h = 0;		//size of the neighborhood considering at first
	unsigned int i = 0;

	srand((unsigned) time(NULL));

//____________________________________________________________________________
// Test1 - Search for OOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (=d)

	t = 3;  
	n = 6;		
	l = 4;		
	d = 1;		
	is_exact = true;	
	m = t*(n-l+1);

	String<DnaString> dataset1;
	appendValue(dataset1,DnaString("ACAGCA"));
	appendValue(dataset1,DnaString("AGGCAG"));
	appendValue(dataset1,DnaString("TCAGTC"));

	//Application of PROJECTION-OOPS
    MotifFinder<Dna, Projection> motif_finder1(t,l,m,d,is_exact);
	findMotif(motif_finder1, dataset1, OOPS());
	//check whether found motif is really an OOPS motif
	SEQAN_TASSERT(isOOPSMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder1,0),
				  d,
				  is_exact)==true);

//____________________________________________________________________________
//
	//Application of ePatternBranching-OOPS
    MotifFinder<Dna, EPatternBranching> motif_finder2(t,l,d,is_exact,h);
	findMotif(motif_finder2, dataset1, OOPS());
	//check whether found motif is really an OOPS motif
	for(i=0; i<length(motif_finder2.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(isOOPSMotif(begin(dataset1),
					  end(dataset1),
				      motif_finder2.set_of_motifs[i],
				      d,
				      is_exact)==true);
	}

//____________________________________________________________________________
// Test2 - Search for OMOPS motifs on a small set of nucleotide sequences
//         given the exact Hamming distance (<=d)

	is_exact = false;

	//Application of PROJECTION-OMOPS
    MotifFinder<Dna, Projection> motif_finder3(t,l,m,d,is_exact);
	findMotif(motif_finder3, dataset1, OMOPS());
	//check whether found motif is really an OMOPS motif
	SEQAN_TASSERT(isOMOPSMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder3, 0),
				  d,
				  
				  is_exact)==true);

//____________________________________________________________________________
//
	//Application of ePatternBranching-OMOPS
	MotifFinder<Dna, EPatternBranching> motif_finder4(t,l,d,is_exact,h);
	findMotif(motif_finder4, dataset1, OMOPS());
	//check whether found motif is really an OMOPS motif
	for(i=0; i<length(motif_finder4.set_of_motifs); ++i)
	{
		SEQAN_TASSERT(isOMOPSMotif(begin(dataset1),
					  end(dataset1),
				      motif_finder4.set_of_motifs[i],
				      d,
				      is_exact)==true);
	}

//____________________________________________________________________________
// Test3 - Search for ZOOPS motifs on a set of small nucleotide sequences
//         given the inexact Hamming distance (<=d)

	//Application of PROJECTION-ZOOPS
    MotifFinder<Dna, Projection> motif_finder5(t,l,m,d,is_exact);
	findMotif(motif_finder5, dataset1, ZOOPS());
	//check whether found motif is really a ZOOPS motif
	SEQAN_TASSERT(isZOOPSMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder5, 0),
				  d,
				  is_exact)==true);

//____________________________________________________________________________
// Test4 - Search for TCM motifs on a set of small nucleotide sequences
//         given the exact Hamming distance (=d)

	is_exact = true;

	//Application of PROJECTION-TCM
    MotifFinder<Dna, Projection> motif_finder6(t,l,m,d,is_exact);
	findMotif(motif_finder6, dataset1, TCM());
	//check whether found motif is really a TCM motif
	SEQAN_TASSERT(isTCMMotif(begin(dataset1),
				  end(dataset1),
				  getMotif(motif_finder6),
				  d,
				  is_exact)==true);
}

int main() 
{
	SEQAN_TREPORT("TEST FIND MOTIF BEGIN")

	Test_exactAlgorithms();
	Test_approximationAlgorithms();

	SEQAN_TREPORT("TEST FIND MOTIF END");

	return 0;
}
