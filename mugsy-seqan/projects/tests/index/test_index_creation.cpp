#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_PROFILE
//#define SEQAN_DEBUG
//#define SEQAN_DEBUG_INDEX

//#define SEQAN_TEST
//#define SEQAN_TEST_SKEW3
//#define SEQAN_TEST_SKEW7

#include "test_index_creation.h"
#include <seqan/index.h>


using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////


bool testIndexCreation()
{
		typedef String<char> TText;
		typedef String<unsigned> TArray;

		TText	text;
		TArray	sa;
		TArray	lcp;
		TArray	child, childExt;
		TText	bwt;

		const int runs = 10;					// conduct 10 test runs 
		const int maxSize = 20 * 1024 * 1024;	// max text size is 20 megabyte
		bool result = true;

		_proFloat timeDelta[12];
		_proFloat timeSum[12];
		for(int i = 0; i < 10; ++i)
			timeSum[i] = 0;
		__int64 textSum = 0;

		static const char* algNames[] = {
			"Skew3        ", 
			"Skew7        ", 
			"ManberMyers  ", 
			"LarssonSadake", 
			"SAQSort      ", 
			"SAQSortQGSR  ", 
			"Skew3Ext     ", 
			"Skew7Ext     ",
			"Kasai        ",
			"KasaiInPlace ",
			"KasaiExt     ",
			"ChildTab     ",
			"ChildTabExt  "
		};

		int TI;
		for(int i = 0; i < runs; ++i) {

			cout << "*** RUN " << i << " ***";
			
			int size = rand() % maxSize;
			TI = 0;

//___randomize_text___________________________________________________________

			resize(text,size);
/*			if (i < runs/2)
				randomize(text);
			else
*/				textRandomize(text);
/*			String<char,External<> > errorText;	// read in text causing an error
			open(errorText,"error.txt");
			text = errorText;
*/
/*			text = "MISSISSIPPI";
			size = length(text);
			cout << "text created (n=" << size << ")" << endl;
*/
			cout << "   textSize: " << length(text) << endl;

//___create_suffix_array______________________________________________________

			resize(sa, size);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, Skew3());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (internal Skew3) failed" << endl;
				result = false;
			}
*/			cout << "."; cout.flush();

			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, Skew7());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (internal Skew7) failed" << endl;
				result = false;
			}
*/			cout << "."; cout.flush();

			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, ManberMyers());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (internal ManberMyers) failed" << endl;
				result = false;
			}
*/			cout << "."; cout.flush();

			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArrayExt(sa, text, LarssonSadakane());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (external LarssonSadakane) failed" << endl;
				result = false;
			}
*/			cout << "."; cout.flush();

			blank(sa);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, SAQSort());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (internal SAQSort) failed" << endl;
				result = false;
			}
			cout << "."; cout.flush();
/*
			blank(sa);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArray(sa, text, QSQGSR(), 3);
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (internal QSQGSR) failed" << endl;
				result = false;
			}
			cout << "."; cout.flush();
*/
			blank(sa);
/*			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArrayExt(sa, text, Skew3());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (external Skew3) failed" << endl;
				result = false;
			}
			cout << "."; cout.flush();

			blank(sa);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createSuffixArrayExt(sa, text, Skew7());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isSuffixArray(sa, text)) {
				cout << "suffix array creation (external Skew7) failed" << endl;
				result = false;
			}
			cout << "."; cout.flush();

//___create_lcp_table_________________________________________________________

			resize(lcp, size);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createLCPTable(lcp, text, sa, KasaiOriginal());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isLCPTable(lcp, sa, text)) {
				cout << "suffix array creation (internal Kasai) failed" << endl;
				result = false;
			}
			cout << "."; cout.flush();

			blank(lcp);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createLCPTable(lcp, text, sa, Kasai());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isLCPTable(lcp, sa, text)) {
				cout << "suffix array creation (internal in-place Kasai) failed" << endl;
				result = false;
			}
			cout << "."; cout.flush();

			blank(lcp);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createLCPTableExt(lcp, text, sa, Kasai());
			timeDelta[TI++] += SEQAN_PROGETTIME;
			if (!isLCPTable(lcp, sa, text)) {
				cout << "suffix array creation (external Kasai) failed" << endl;
				result = false;
			}
			cout << "."; cout.flush();

//___create_child_table_______________________________________________________

			resize(child, size);
			for(int i=0; i<size; ++i)
				child[i] = supremumValue<unsigned>();
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createChildTable(child, lcp);
			timeDelta[TI++] += SEQAN_PROGETTIME;
			cout << "."; cout.flush();

			unsigned undefs=0;
			for(int i=0; i<size; ++i)
				if (child[i] == supremumValue<unsigned>()) ++undefs;
			if (undefs) ::std::cout << undefs << " undefined values";

			resize(childExt, size);
			timeDelta[TI] = -SEQAN_PROGETTIME;
			createChildTableExt(childExt, lcp);
			timeDelta[TI++] += SEQAN_PROGETTIME;
			cout << "."; cout.flush();

			if (!isEqual(child, childExt)) {
				cout << "child table creation failed" << endl;
				result = false;
			}
*/
//___update_performance_table_________________________________________________

			for(int i=0; i<TI; ++i) {
				timeSum[i] += timeDelta[i];
				textSum += length(text);
			}

			cout << " OK!" << endl;

		}
		cout << "*** TIME RESULTS (sec/MB) ***" << endl;
		for(int i=0; i<TI; ++i)
			cout << algNames[i] << " " << 1024.0*1024.0 * timeSum[i] / textSum << endl;

		return result;
}

/*
int main() {
	return (testIndexCreation())? 0: 1;
}
*/

