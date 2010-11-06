#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST
#define TEST_PATH "projects/tests/index/"

#include <seqan/index.h>
#include <seqan/file.h>

using namespace std;
using namespace seqan;


void testGappedShapes()
{
	String<char> shape_string = "0010011011101";
	Shape<Dna,GenericShape> shape1;
	stringToShape(shape1, shape_string);
	Shape<Dna,GenericShape> shape2 = Shape<Dna,GenericShape>(shape1);

	SEQAN_ASSERT(shape1.weight == shape2.weight);
	SEQAN_ASSERT(shape1.span == shape2.span);
	SEQAN_ASSERT(shape1.diffs == shape2.diffs);
	SEQAN_ASSERT(length(shape1) == length(shape2));
	SEQAN_ASSERT(weight(shape1) == weight(shape2));
/*
	Shape<Dna,GenericShape> shape3 = Shape<Dna,GenericShape>(5, 13);
	shape3[0]=2;
	shape3[1]=2;
	shape3[2]=1;
	shape3[3]=1;
	shape3[4]=2;
	shape3[5]=1;
	shape3[6]=1;
	shape3[7]=2;
	for(int i = 0; i < 8; ++i)
        SEQAN_ASSERT(shape1[i] == shape3[i]);
*/
}


void testUngappedShapes()
{
	Shape<Dna,SimpleShape> shape1;
	resize(shape1, 4);
	Shape<Dna,SimpleShape> shape2 = Shape<Dna,SimpleShape>(shape1);

	SEQAN_ASSERT(shape1.span == shape2.span);
	SEQAN_ASSERT(shape1.leftFactor == shape2.leftFactor);
	SEQAN_ASSERT(length(shape1) == length(shape2));
	SEQAN_ASSERT(weight(shape1) == weight(shape2));
	

	Shape<Dna,SimpleShape> shape3 = Shape<Dna,SimpleShape>(4);
	SEQAN_ASSERT(shape3.leftFactor == 64);
	SEQAN_ASSERT(shape1.leftFactor == shape3.leftFactor);


}

/*
void testQGramIndexSchnell()
{
	clock_t start, finish;
	double duration;

	String<Dna> text;
	fstream strm_t;
	strm_t.open(TEST_PATH "fasta.txt", ios_base::in);
	read(strm_t, text, Fasta());
	String<Dna> next;
	resize(next,length(text));
	for(int i = 0; i < 1; ++i)							// datei zu ende?
	{
		arrayCopyForward(begin(text),end(text),begin(next));
		append(text, next);
	}
	strm_t.close();

	String<char> shape_string = "1000100100001110011";
	int q = length(shape_string);
	Shape<Dna,GenericShape> shape;
	stringToShape(shape, shape_string);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);	
	
	String<TPosition> pos;	
	int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(shape));
	pos_size += 1;	
	resize(pos, pos_size);

	start = clock();
	createQGramIndex(index, pos, text, shape);
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "\nQGramIndex bauen dauert: " << duration << " Sekunden.\n\n";
	
	
}
*/
/*
void testGappedQGramIndex()
{
	String<Dna> text = "CTGAACCCTAAACCCT";
	String<char> shape_string = "101";
	int q = length(shape_string);
	Shape<Dna,GenericShape> shape;
	stringToShape(shape, shape_string);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 2);
	
	String<TPosition> pos;
    int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, weight(shape));
	pos_size += 1;	
	resize(pos, pos_size);

	createQGramIndex(index, pos, text, shape);
	
	SEQAN_ASSERT(pos[0] == 0);
	SEQAN_ASSERT(pos[1] == 1);
	SEQAN_ASSERT(pos[2] == 5);
	SEQAN_ASSERT(pos[3] == 5);
	SEQAN_ASSERT(pos[4] == 5);
	SEQAN_ASSERT(pos[5] == 6);
	SEQAN_ASSERT(pos[6] == 8);
	SEQAN_ASSERT(pos[7] == 9);
	SEQAN_ASSERT(pos[8] == 11);
	SEQAN_ASSERT(pos[9] == 12);
	SEQAN_ASSERT(pos[10] == 12);
	SEQAN_ASSERT(pos[11] == 12);
	SEQAN_ASSERT(pos[12] == 12);
	SEQAN_ASSERT(pos[13] == 14);
	SEQAN_ASSERT(pos[14] == 14);
	SEQAN_ASSERT(pos[15] == 14);
	SEQAN_ASSERT(pos[16] == 14);

	SEQAN_ASSERT(index[0] == 9);
	SEQAN_ASSERT(index[1] == 11);
	SEQAN_ASSERT(index[2] == 10);
	SEQAN_ASSERT(index[3] == 4);
	SEQAN_ASSERT(index[4] == 3);
	SEQAN_ASSERT(index[5] == 7);
	SEQAN_ASSERT(index[6] == 12);
	SEQAN_ASSERT(index[7] == 5);
	SEQAN_ASSERT(index[8] == 0);
	SEQAN_ASSERT(index[9] == 13);
	SEQAN_ASSERT(index[10] == 6);
	SEQAN_ASSERT(index[11] == 2);
	SEQAN_ASSERT(index[12] == 8);
	SEQAN_ASSERT(index[13] == 1);
	
}
*/
void testUngappedQGramIndex()
{
	String<Dna> text = "CTGAACCCTAAACCCT";
	int q = 2;
	Shape<Dna,SimpleShape> shape;
	resize(shape, q);

	typedef Position<String<Dna> >::Type TPosition;
	String<TPosition> index;
	resize(index, length(text) - q + 1);
	
	String<TPosition> pos;
    int pos_size = _intPow((unsigned)ValueSize<Dna>::VALUE, q);
	pos_size += 1;	
	resize(pos, pos_size);

	createQGramIndex(index, pos, text, shape);
	
	
	SEQAN_ASSERT(pos[0] == 0);
	SEQAN_ASSERT(pos[1] == 3);
	SEQAN_ASSERT(pos[2] == 5);
	SEQAN_ASSERT(pos[3] == 5);
	SEQAN_ASSERT(pos[4] == 5);
	SEQAN_ASSERT(pos[5] == 5);
	SEQAN_ASSERT(pos[6] == 9);
	SEQAN_ASSERT(pos[7] == 9);
	SEQAN_ASSERT(pos[8] == 12);
	SEQAN_ASSERT(pos[9] == 13);
	SEQAN_ASSERT(pos[10] == 13);
	SEQAN_ASSERT(pos[11] == 13);
	SEQAN_ASSERT(pos[12] == 13);
	SEQAN_ASSERT(pos[13] == 14);
	SEQAN_ASSERT(pos[14] == 14);
	SEQAN_ASSERT(pos[15] == 15);

	SEQAN_ASSERT(index[0] == 3);
	SEQAN_ASSERT(index[1] == 9);
	SEQAN_ASSERT(index[2] == 10);
	SEQAN_ASSERT(index[3] == 4);
	SEQAN_ASSERT(index[4] == 11);
	SEQAN_ASSERT(index[5] == 5);
	SEQAN_ASSERT(index[6] == 6);
	SEQAN_ASSERT(index[7] == 12);
	SEQAN_ASSERT(index[8] == 13);
	SEQAN_ASSERT(index[9] == 0);
	SEQAN_ASSERT(index[10] == 7);
	SEQAN_ASSERT(index[11] == 14);
	SEQAN_ASSERT(index[12] == 2);
	SEQAN_ASSERT(index[13] == 8);
	SEQAN_ASSERT(index[14] == 1);
}

//////////////////////////////////////////////////////////////////////////////

void testQGramFind()
{
	typedef Index<String<char>, Index_QGram<UngappedShape<2> > > TQGramIndex;
	TQGramIndex idx("to be or not to be");
	Finder<TQGramIndex> finder(idx);

	SEQAN_TASSERT(find(finder, "be"))
	SEQAN_TASSERT(position(finder) == 3);
	SEQAN_TASSERT(find(finder, "be"))
	SEQAN_TASSERT(position(finder) == 16);
	SEQAN_TASSERT(!find(finder, "be"))
/*
	while (find(finder, "be"))
	{
		cout << position(finder) << "\n";
	}
*/
}

//////////////////////////////////////////////////////////////////////////////

void testQGramIndex()
{
	SEQAN_TREPORT("TEST QGRAM BEGIN")

//	testGappedQGramIndex();
	testUngappedQGramIndex();
//	testQGramIndexSchnell();
	testGappedShapes();
	testUngappedShapes();

	testQGramFind();

	debug::verifyCheckpoints("projects/library/seqan/index/index_qgram.h");
	debug::verifyCheckpoints("projects/library/seqan/index/shape_base.h");
	debug::verifyCheckpoints("projects/library/seqan/index/shape_gapped.h");

	debug::verifyCheckpoints("projects/library/seqan/index/index_qgram_find.h");


	SEQAN_TREPORT("TEST QGRAM END")
}
