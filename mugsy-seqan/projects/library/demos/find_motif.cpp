///A tutorial about motif search.
#include <iostream>
#include "seqan/find_motif.h"

using namespace seqan;

///Function to output found motifs.
template <typename TMotifFinder>
void printMotifs(TMotifFinder & finder)
{
	for (int i = 0; i < (int) motifCount(finder); ++i)
	{
		::std::cout << i << ": " << getMotif(finder, i) << ::std::endl;
	}
}

int main() 
{
	::std::srand((unsigned) time(NULL));

///Motif search on a small set of nucleotide sequences.
	unsigned int t = 3;		//number of input sequences
	unsigned int n = 6;		//length of sequence
	unsigned int l = 4;		//length of motif
	unsigned int d = 1;		//number of substitutions
	bool is_exact = true;	//size of Hamming distance
	unsigned int h = 0;		//size of the neighborhood considering at first

	String<DnaString> dataset;
	appendValue(dataset,DnaString("ACAGCA"));
	appendValue(dataset,DnaString("AGGCAG"));
	appendValue(dataset,DnaString("TCAGTC"));
	
///Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb1(t,l,d,is_exact,h);
	findMotif(finder_epb1,dataset,OMOPS());
	::std::cout << getMotif(finder_epb1) << ::std::endl;

///Application of ePatternBranching (h=0)
	MotifFinder<Dna, EPatternBranching> finder_epb2(t,l,d,is_exact,h);
	findMotif(finder_epb2,dataset,OOPS());
	::std::cout << getMotif(finder_epb2) << ::std::endl;

///Application of PMS1-ZOOPS 
	MotifFinder<Dna, PMS1> finder_pms1(l,d,is_exact);
	findMotif(finder_pms1,dataset,ZOOPS());
	printMotifs(finder_pms1); 

///Application of PMSP-TCM
	MotifFinder<Dna, PMSP> finder_pmsp(l,d,is_exact);
	findMotif(finder_pmsp,dataset,TCM());
	printMotifs(finder_pmsp); 
	
///Application of PROJECTION-OOPS
	unsigned int m = t*(n-l+1);
    MotifFinder<Dna, Projection> finder_proj(t,l,m,d,is_exact);
	findMotif(finder_proj, dataset, OOPS());
	printMotifs(finder_proj);

///Application of PROJECTION-OMOPS
    MotifFinder<Dna, Projection> finder_proj_omops(t,l,m,d,is_exact);
	findMotif(finder_proj_omops, dataset, OMOPS());
	printMotifs(finder_proj_omops);

///Application of PROJECTION-ZOOPS
	MotifFinder<Dna, Projection> finder_proj_zoops(t,l,m,d,is_exact);
	findMotif(finder_proj_zoops, dataset, ZOOPS());
	printMotifs(finder_proj_zoops);
	
///Application of PROJECTION-TCM
    MotifFinder<Dna, Projection> finder_proj_tcm(t,l,m,d,is_exact);
	findMotif(finder_proj_tcm, dataset, TCM());
	printMotifs(finder_proj_tcm);

	return 0;
}

