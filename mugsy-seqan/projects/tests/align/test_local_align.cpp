#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/align.h>


using namespace std;
using namespace seqan;


void testLocalAlign(){

	//align two sequences using Smith-Waterman-algorithm
	String<char> str0 = "ataagcgtctcg";
	String<char> str1 = "tcatagagttgc";
	
	Align< String<char>, ArrayGaps> ali;
	resize(rows(ali), 2);
	setSource(row(ali, 0), str0);
	setSource(row(ali, 1), str1);

	Score<int> score_type = Score<int>(2,-1,-2,0) ;
	LocalAlignmentFinder<int> sw_finder = LocalAlignmentFinder<int>(ali);
	
	int cutoff = 0;

	int score = smithWaterman(ali,sw_finder,score_type,cutoff);
	
	SEQAN_TASSERT(score == 9);
	SEQAN_TASSERT(row(ali,0) == "ataagcgt");
	SEQAN_TASSERT(row(ali,1) == "ata-gagt");
	
	int i = 1;
	while (true){
		
		score = smithWatermanGetNext(ali,sw_finder,score_type,cutoff);
		if(score==0){
	//		cout <<"No more alignments satisfying score > "<<cutoff<<"found.\n";		
			break;
		}
		if(i == 1){
			SEQAN_TASSERT(score == 5);
			SEQAN_TASSERT(row(ali,0) == "tc-tcg");
			SEQAN_TASSERT(row(ali,1) == "tcatag");
		}
		if(i == 2){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "taagcgtctcg");
			SEQAN_TASSERT(row(ali,1) == "tcatagagttg");
		}
		if(i == 3){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "gc");
			SEQAN_TASSERT(row(ali,1) == "gc");
		}
		if(i == 4){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "ag");
			SEQAN_TASSERT(row(ali,1) == "ag");
		}
		if(i == 5){
			SEQAN_TASSERT(score == 4);
			SEQAN_TASSERT(row(ali,0) == "tc");
			SEQAN_TASSERT(row(ali,1) == "tc");
		}
		if(i == 6){
			SEQAN_TASSERT(score == 3);
			SEQAN_TASSERT(row(ali,0) == "cgt");
			SEQAN_TASSERT(row(ali,1) == "cat");
		}
		if(i == 7){
			SEQAN_TASSERT(score == 3);
			SEQAN_TASSERT(row(ali,0) == "ata");
			SEQAN_TASSERT(row(ali,1) == "aga");
		}
		if(i == 8){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "g");
			SEQAN_TASSERT(row(ali,1) == "g");
		}
		if(i == 9){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		if(i == 10){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "a");
			SEQAN_TASSERT(row(ali,1) == "a");
		}
		if(i == 11){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "c");
			SEQAN_TASSERT(row(ali,1) == "c");
		}
		if(i == 12){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "c");
			SEQAN_TASSERT(row(ali,1) == "c");
		}
		if(i == 13){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "a");
			SEQAN_TASSERT(row(ali,1) == "a");
		}
		if(i == 14){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "a");
			SEQAN_TASSERT(row(ali,1) == "a");
		}
		if(i == 15){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		if(i == 16){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "g");
			SEQAN_TASSERT(row(ali,1) == "g");
		}
		if(i == 17){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		if(i == 18){
			SEQAN_TASSERT(score == 2);
			SEQAN_TASSERT(row(ali,0) == "t");
			SEQAN_TASSERT(row(ali,1) == "t");
		}
		++i;
	}
	
//test if every cell has been reduced to 0 
//only makes sense if cutoff=0
	if(cutoff==0){
		int str0len = length(str0) + 1;
		int str1len = length(str1) + 1;
		bool check = true;
		for(int i = 0; i <str1len; ++i){
			for(int j=0;j<str0len;++j){
				if(getValue(sw_finder.matrix_,(i*str0len)+j)!=0){
					check = false;
				}
			}
		}

		SEQAN_TASSERT(check == true);

	}

//desweiteren nur so:
	push(sw_finder.pq_,LocalAlignmentFinder<int>::TPQEntry());
	SEQAN_TASSERT(empty(sw_finder.pq_) == false);
	clear(sw_finder.pq_);
	SEQAN_TASSERT(empty(sw_finder.pq_) == true);




}


void testLocalAlign2()
{

//new interface

	String<char> str0 = "ataagcgtctcg";
	String<char> str1 = "tcatagagttgc";
	
	Align< String<char>, ArrayGaps> ali;
	resize(rows(ali), 2);
	setSource(row(ali, 0), str0);
	setSource(row(ali, 1), str1);

	Score<int> score_type = Score<int>(2,-1,-2,0) ;
	LocalAlignmentFinder<int> sw_finder = LocalAlignmentFinder<int>(ali);
	
	int score = localAlignment(ali, sw_finder, score_type, 4);
	SEQAN_TASSERT(score == 9);
	SEQAN_TASSERT(row(ali,0) == "ataagcgt");
	SEQAN_TASSERT(row(ali,1) == "ata-gagt");
	
	score = localAlignment(ali, sw_finder, score_type, 4);
	SEQAN_TASSERT(score == 5);
	SEQAN_TASSERT(row(ali,0) == "tc-tcg");
	SEQAN_TASSERT(row(ali,1) == "tcatag");

	score = localAlignment(ali, sw_finder, score_type, 4, WatermanEggert());
	SEQAN_TASSERT(score == 0);
}



void Main_TestLocalAlign() 
{

	SEQAN_TREPORT("TEST LOCAL ALIGN BEGIN")

	testLocalAlign();
	testLocalAlign2();



	debug::verifyCheckpoints("projects/library/seqan/misc/priority_type_base.h");
	debug::verifyCheckpoints("projects/library/seqan/misc/priority_type_heap.h");
	debug::verifyCheckpoints("projects/library/seqan/align/align_local_dynprog.h");

		 
	SEQAN_TREPORT("TEST LOCAL ALIGN END")
}
