#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/file.h>
#include <seqan/align.h>

#include "test_embl.h"


#define TEST_PATH "projects/tests/file/"
#define LIB_PATH "projects/library/seqan/file/"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////


void Test_Stream()
{
	//string => stream
	String<char> str_1("The first test file string.");

	fstream strm_1;
	strm_1.open(TEST_PATH "testfile.txt", ios_base::out | ios_base::trunc);

	strm_1 << str_1;
	strm_1.close();

	//assign stream => string
	String<char> str_3 = "";
	fstream strm_3;
	strm_3.open(TEST_PATH "testfile.txt", ios_base::in);
	strm_3 >> str_3;
	strm_3.close();
	SEQAN_TASSERT(str_3 == str_1)


	//const string => stream
	fstream strm_4;
	strm_4.open(TEST_PATH "testfile.txt", ios_base::out | ios_base::trunc);
	strm_4 << "Testfile";
	strm_4.close();

	String<char> str_5 = "";
	fstream strm_5;
	strm_5.open(TEST_PATH "testfile.txt", ios_base::in);
	strm_5 >> str_5;
	strm_5.close();
	SEQAN_TASSERT(str_5 == "Testfile")


//	String<char> str_4("This test file was made by " __FILE__ " too.\n");

/*
	str_2 = "";
	fstream strm_3;
	strm_3.open(TEST_PATH "testfile.txt", ios_base::in);
	assign(str_2, strm_3, Insist());
	strm_3.close();
	SEQAN_TASSERT(str_2 == str_1)

	str_2 = "";
	fstream strm_4;
	strm_4.open(TEST_PATH "testfile.txt", ios_base::in);
	assign(str_2, strm_4, 9, Insist());
	strm_4.close();
	SEQAN_TASSERT(str_2 == "This test")

	String<char> str_3;
	fstream strm_5;
	strm_5.open(TEST_PATH "testfile.txt", ios_base::in);
	assign(str_3, strm_5, 9, Generous());
	strm_5.close();
	SEQAN_TASSERT(str_3 == "This test")

	str_2 = "";
	fstream strm_6;
	strm_6.open(TEST_PATH "testfile.txt", ios_base::in);
	assign(str_2, strm_6, Limit());
	strm_6.close();
	SEQAN_TASSERT(str_2 == str_1)

	str_2 = "";
	fstream strm_7;
	strm_7.open(TEST_PATH "testfile.txt", ios_base::in);
	assign(str_2, strm_7, 9, Limit());
	strm_7.close();
	SEQAN_TASSERT(str_2 == "This test")

	char str_4[200];
	fstream strm_8;
	strm_8.open(TEST_PATH "testfile.txt", ios_base::in);
	assign(str_4, strm_8);
	strm_8.close();
	SEQAN_TASSERT(str_4 == str_1)

	fstream strm_9;
	strm_9.open(TEST_PATH "testfile.txt", ios_base::in);
	assign(str_4, strm_9, 9);
	strm_9.close();
	SEQAN_TASSERT(isEqual(str_4, "This test"))

//____________________________________________________________________________
//append stream => string
	String<char> str_5("Start");
	append(str_5, str_1);

	String<char> str_6("Start");
	fstream strm_10;
	strm_10.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_6, strm_10);
	strm_10.close();
	SEQAN_TASSERT(str_6 == str_5)

	String<char> str_7("Start");
	fstream strm_11;
	strm_11.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_7, strm_11, 14, Generous());
	strm_11.close();
	SEQAN_TASSERT(str_7 == "StartThis test")

	String<char> str_8("Start");
	fstream strm_11a;
	strm_11a.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_8, strm_11a, 200, Generous());
	strm_11a.close();
	SEQAN_TASSERT(str_8 == str_5)

	str_2 = "Start";
	fstream strm_12;
	strm_12.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_2, strm_12, Insist());
	strm_12.close();
	SEQAN_TASSERT(str_2 == str_5)

	str_2 = "Start";
	fstream strm_13;
	strm_13.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_2, strm_13, 14, Insist());
	strm_13.close();
	SEQAN_TASSERT(str_2 == "StartThis test")

	str_2 = "Start";
	fstream strm_14;
	strm_14.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_2, strm_14, Limit());
	strm_14.close();
	SEQAN_TASSERT(str_2 == str_5)

	str_2 = "Start";
	fstream strm_15;
	strm_15.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_2, strm_15, 14, Limit());
	strm_15.close();
	SEQAN_TASSERT(str_2 == "StartThis test")

	assign(str_4, "Start");
	fstream strm_16;
	strm_16.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_4, strm_16);
	strm_16.close();
	SEQAN_TASSERT(str_4 == str_5)

	assign(str_4, "Start");
	fstream strm_17;
	strm_17.open(TEST_PATH "testfile.txt", ios_base::in);
	append(str_4, strm_17, 14);
	strm_17.close();
	SEQAN_TASSERT(isEqual(str_4, "StartThis test"))

//____________________________________________________________________________
//replace stream => string
	String<char> str_9("start ");
	append(str_9, str_1);
	append(str_9, " end");

	String<char> str_10("start middle end");
	fstream strm_18;
	strm_18.open(TEST_PATH "testfile.txt", ios_base::in);
	replace(str_10, 6, 12, strm_18);
	strm_18.close();
	SEQAN_TASSERT(str_10 == str_9)

	str_10 = "start middle end";
	fstream strm_19;
	strm_19.open(TEST_PATH "testfile.txt", ios_base::in);
	replace(str_10, 6, 12, strm_19, 15);
	strm_19.close();
	SEQAN_TASSERT(str_10 == "start This test")

	str_10 = "start middle end";
	fstream strm_20;
	strm_20.open(TEST_PATH "testfile.txt", ios_base::in);
	replace(str_10, 6, 12, strm_20, Limit());
	strm_20.close();
	SEQAN_TASSERT(str_10 == str_9)

	str_10 = "start middle end";
	fstream strm_21;
	strm_21.open(TEST_PATH "testfile.txt", ios_base::in);
	replace(str_10, 6, 12, strm_21, 15, Limit());
	strm_21.close();
	SEQAN_TASSERT(str_10 == "start This test")

	assign(str_4, "start middle end");
	fstream strm_22;
	strm_22.open(TEST_PATH "testfile.txt", ios_base::in);
	replace(str_4, 6, 12, strm_22);
	strm_22.close();
	SEQAN_TASSERT(str_4 == str_9)

	assign(str_4, "start middle end");
	fstream strm_23;
	strm_23.open(TEST_PATH "testfile.txt", ios_base::in);
	replace(str_4, 6, 12, strm_23, 15);
	strm_23.close();
	SEQAN_TASSERT(isEqual(str_4, "start This test"))

//____________________________________________________________________________
//Test limited string => stream

	fstream strm_24;
	strm_24.open(TEST_PATH "testfile.txt", ios_base::out | ios_base::trunc);
	write(strm_24, str_1, 9);
	strm_24.close();

	fstream strm_25;
	strm_25.open(TEST_PATH "testfile.txt", ios_base::in);
	strm_25 >> str_2;
	strm_25.close();
	SEQAN_TASSERT(isEqual(str_2, "This test"))

//____________________________________________________________________________
//Test stream => segment

	str_2 = "begin middle end";
	Segment<String<char> > seg_1(str_2, 6, 12);
	fstream strm_26;
	strm_26.open(TEST_PATH "testfile.txt", ios_base::in);
	strm_26 >> seg_1;
	strm_26.close();
	SEQAN_TASSERT(isEqual(str_2, "begin This test end"))

	str_2 = "begin middle end";
	Segment<String<char> > seg_2(str_2, 6, 12);
	fstream strm_27;
	strm_27.open(TEST_PATH "testfile.txt", ios_base::in);
	append(seg_2, strm_27);
	strm_27.close();
	SEQAN_TASSERT(isEqual(str_2, "begin middleThis test end"))

	str_2 = "begin one two three end";
	Segment<String<char> > seg_3(str_2, 6, 19);
	fstream strm_28;
	strm_28.open(TEST_PATH "testfile.txt", ios_base::in);
	replace(seg_3, 4, 7, strm_28);
	strm_28.close();
	SEQAN_TASSERT(isEqual(str_2, "begin one This test three end"))

//____________________________________________________________________________
*/
}

//////////////////////////////////////////////////////////////////////////////

void Test_CStream()
{
//____________________________________________________________________________
//string => stream: create testfile

	String<char> str_1("This test file was made by " __FILE__ ".\n");

	FILE * file_1 = fopen(TEST_PATH "testfile.txt", "w");
	file_1 << str_1;
	fclose(file_1);

/*
//assign FILE => string
//____________________________________________________________________________

	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	String<char> str_2;
	assign(str_2, file_1);
	fclose(file_1);
	SEQAN_TASSERT(str_2 == str_1)

	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	assign(str_2, file_1, 9);
	fclose(file_1);
	SEQAN_TASSERT(str_2 == "This test")

	char str_3[200]="";
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	assign(str_3, file_1);
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_3, str_1))

	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	assign(str_3, file_1, 9);
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_3, "This test"))

//____________________________________________________________________________
//append FILE => string

	String<char> str_5("Start");
	append(str_5, str_1);

	assign(str_2, "Start");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	append(str_2, file_1);
	fclose(file_1);
	SEQAN_TASSERT(str_2 == str_5)

	assign(str_2, "Start");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	append(str_2, file_1, 14);
	fclose(file_1);
	SEQAN_TASSERT(str_2 == "StartThis test")

	assign(str_3, "Start");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	append(str_3, file_1);
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_3, str_5))

	assign(str_3, "Start");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	append(str_3, file_1, 14);
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_3, "StartThis test"))

//____________________________________________________________________________
//replace FILE => string

	String<char> str_9("start ");
	append(str_9, str_1);
	append(str_9, " end");

	assign(str_2, "start middle end");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	replace(str_2, 6, 12, file_1);
	fclose(file_1);
	SEQAN_TASSERT(str_2 == str_9)

	assign(str_2, "start middle end");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	replace(str_2, 6, 12, file_1, 15);
	fclose(file_1);
	SEQAN_TASSERT(str_2 == "start This test")

	assign(str_3, "start middle end");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	replace(str_3, 6, 12, file_1);
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_3, str_9))

	assign(str_3, "start middle end");
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	replace(str_3, 6, 12, file_1, 15);
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_3, "start This test"))
	
//____________________________________________________________________________
// assign string => FILE

	str_1 = "first seqan string test";
	FILE * file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, str_1, 18);
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, infix(str_1, 0, 18)))

	assign(str_3, "first char array test");
	FILE * file_2 = fopen("testfile2.txt", "w");
	write(file_2, str_3);
	fclose(file_2);
	file_1 = fopen("testfile2.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, str_3))

	assign(str_3, "second char array test");
	file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, str_3, 17);
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, infix(str_3, 0, 17)))

	char const * str_20 = "first const char array test";
	file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, str_20);
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, str_20))

	char const * str_21 = "second const char array test";
	file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, str_21, 23);
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, infix(str_21, 0, 23)))

	String<char> str_22("begin middle end");
	Segment<String<char> > infix_1(str_22, 6, 12);
	file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, infix_1);
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, infix_1))

	assign(infix_1, "prefix suffix");
	file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, infix_1, 6);
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, infix(infix_1, 0, 6)))

	file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, 'a');
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, 'a'))

	file_2 = fopen(TEST_PATH "testfile.txt", "w");
	write(file_2, 'a', 6);
	fclose(file_2);
	file_1 = fopen(TEST_PATH "testfile.txt", "r");
	file_1 >> str_2;
	fclose(file_1);
	SEQAN_TASSERT(isEqual(str_2, 'a'))
*/

//____________________________________________________________________________

}
//////////////////////////////////////////////////////////////////////////////

void Test_Raw()
{
//____________________________________________________________________________

	//write
	FILE * file_3 = fopen(TEST_PATH "testfile.txt", "w");
	String<char> str_3("this is a test string"); 
	write(file_3, str_3);
	fclose(file_3);

	FILE * file_4 = fopen(TEST_PATH "testfile.txt", "r");
	SEQAN_TASSERT(file_4);

	//readID
	String<char> str_4("init");
	readID(file_4, str_4, Raw() );
	SEQAN_TASSERT(str_4 == "");

	//read
	read(file_4, str_4);
	SEQAN_TASSERT(str_4 == str_3);
	fclose(file_4);

	//read limit
	FILE * file_4b = fopen(TEST_PATH "testfile.txt", "r");
	SEQAN_TASSERT(file_4b);

	read(file_4b, str_4, 4);
	SEQAN_TASSERT(str_4 == "this");
	fclose(file_4b);

//____________________________________________________________________________
// raw and stingstream

	str_3 = "hello you fellow";
	stringstream ss_1;
	write(ss_1, str_3);

	read(ss_1, str_4, Raw());
	SEQAN_TASSERT(str_4 == str_3);


	str_3 = "gogogo dududu";
	stringstream ss_2;
	write(ss_2, str_3);
	read(ss_2, str_4, 5);
	SEQAN_TASSERT(str_4 == "gogog");


//	::std::cout << sstream_1;

/*
	//string => stream
	String<char> str_1("The first test file string.");

	fstream strm_1;
	strm_1.open(TEST_PATH "testfile.txt", ios_base::out | ios_base::trunc);
	strm_1 << str_1;
	strm_1.close();

	//assign stream => string
	String<char> str_3 = "";
	fstream strm_3;
	strm_3.open(TEST_PATH "testfile.txt", ios_base::in);
	strm_3 >> str_3;
	strm_3.close();
	SEQAN_TASSERT(str_3 == str_1)
*/

//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

void Test_Fasta_Read(char * path)
{
//____________________________________________________________________________
// FASTA from C stream

	FILE * file_1 = fopen(path, "rb");
	SEQAN_TASSERT(file_1);

	String<Dna> str_1;
	read(file_1, str_1, Fasta());
	SEQAN_TASSERT(str_1 == "ACGT");

	String<char> str_2;

	read(file_1, str_2, Fasta());
	SEQAN_TASSERT(length(str_2) == 0);

	read(file_1, str_2, Fasta());
	SEQAN_TASSERT(length(str_2) == 1065);

	read(file_1, str_2, Fasta());
	SEQAN_TASSERT(length(str_2) == 200);

	fclose(file_1);

//____________________________________________________________________________
// FASTA from iostream

	fstream strm_1;
	strm_1.open(path, ios_base::in | ios_base::binary);
	read(strm_1, str_2, Fasta());
	strm_1.close();
	SEQAN_TASSERT(str_2 == "ACGT");
}

void Test_Fasta_Write()
{
//____________________________________________________________________________
// FASTA to C stream

	FILE * file_3 = fopen(TEST_PATH "my_fasta.txt", "wb");
	String<Dna> str_3("acgtufacgtufacgtufacgtufacgtufacgtufacgtufacgtufaaaaaaaaaauuuuuuuuucccccccccccggggg"); 
	write(file_3, str_3, "Identifier1", Fasta());
	write(file_3, "", "Empty Entry", Fasta());
	write(file_3, str_3, "Identifier2", Fasta());
	fclose(file_3);

	//test it
	FILE * file_4 = fopen(TEST_PATH "my_fasta.txt", "rb");
	SEQAN_TASSERT(file_4);

	String<char> str_4;
	readID(file_4, str_4, Fasta());
	SEQAN_TASSERT(str_4 == "Identifier1");

	String<Dna> str_4a;
	read(file_4, str_4a, Fasta());
	SEQAN_TASSERT(str_4a == str_3);

	String<Dna> str_empty;
	read(file_4, str_empty, Fasta());
	SEQAN_TASSERT(length(str_empty) == 0)

	read(file_4, str_4a, Fasta());
	SEQAN_TASSERT(str_4a == str_3);

//____________________________________________________________________________
// Test virtual file format object //outdated!!
/*
	FileFormat<FILE *, String<Dna>, Fasta> fasta_format;
	FileFormat<FILE *, String<Dna> > & ff = fasta_format;

	//write
	FILE * file_5 = fopen(TEST_PATH "my_fasta.txt", "w");
	String<Dna> str_5("acgt"); 
	write(file_5, str_5, "Identifier3", ff);
	write(file_5, str_5, "Identifier4", ff);
	write(file_5, str_5, "Identifier5", ff);
	write(file_5, str_5, "Identifier6", ff);
	fclose(file_5);

	//readID
	FILE * file_6 = fopen(TEST_PATH "my_fasta.txt", "r");
	SEQAN_TASSERT(file_6);

	String<char> str_6;
	readID(file_6, str_6, ff);
	SEQAN_TASSERT(str_6 == "Identifier3");

	char str_6b[100];
	readID(file_6, str_6b, ff);
	SEQAN_TASSERT(isEqual(str_6b, "Identifier3"));

	//read
	String<Dna> str_7;
	read(file_6, str_7, ff);
	SEQAN_TASSERT(str_7 == str_5);

	read(file_6, str_7, 2, ff);
	SEQAN_TASSERT(str_7 == "ac");

	readID(file_6, str_6, ff);
	SEQAN_TASSERT(str_6 == "Identifier5");

	//goNext
	goNext(file_6, ff);

	readID(file_6, str_6, ff);
	SEQAN_TASSERT(str_6 == "Identifier6");

	//comparison
	SEQAN_TASSERT(ff == fasta_format);
	SEQAN_TASSERT(ff == Fasta());
	SEQAN_TASSERT(Fasta() == ff);

	SEQAN_TASSERT(!(ff != fasta_format));
	SEQAN_TASSERT(ff != Raw());
	SEQAN_TASSERT(Raw() != ff);
*/
}

//////////////////////////////////////////////////////////////////////////////

void Test_FastaAlign() {
//____________________________________________________________________________
// FASTA Align from C stream

	FILE * file_1 = fopen(TEST_PATH "fasta_align.txt", "r");
	SEQAN_TASSERT(file_1);
	
	Align<String<char>, ArrayGaps> align;
	read(file_1, align, FastaAlign());
	SEQAN_TASSERT(source(row(align,0)) == "MKVILLFVLAVFTVFVSSRGIPPEEQSQFLEFQDKFNKKYSHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDEFKNYYLNNKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGPLAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII");
	SEQAN_TASSERT(length(source(row(align,1))) == 362);
	SEQAN_TASSERT(length(source(row(align,2))) == 335);
	String<String<char> > str_ids1;
	readIDs(file_1, str_ids1, FastaAlign());
	SEQAN_TASSERT(value(str_ids1,0) == "CYS1_DICDI");
	SEQAN_TASSERT(value(str_ids1,1) == "ALEU_HORVU");
	SEQAN_TASSERT(value(str_ids1,2) == "CATH_HUMAN");

	fclose(file_1);

	FILE * file_11 = fopen(TEST_PATH "fasta_align_dna.txt", "r");
	SEQAN_TASSERT(file_11);
	
	Align<String<Dna>, ArrayGaps> align11;
	read(file_11, align11, FastaAlign());
	SEQAN_TASSERT(source(row(align11,0)) == "CTACGAAAGGTCGTGTCACGATGTCCGCAAGGGATGGCATTGCATAGAGGAATTGATTGCAACCTACGAAA");
	SEQAN_TASSERT(source(row(align11,1)) == "CTTAATGTCCCGCGTACAAGGGATAGCATGTGGCATAGAGGAATAGAATAGCAGCCTACGAAA");
	String<String<char> > str_ids11;
	readIDs(file_11, str_ids11, FastaAlign());
	SEQAN_TASSERT(value(str_ids11,0) == "SEQ1");
	SEQAN_TASSERT(value(str_ids11,1) == "SEQ2");
	fclose(file_11);

	FILE * file_44 = fopen(TEST_PATH "my_fasta_align_dna.txt", "w");
	write(file_44, align11, str_ids11, FastaAlign());
	fclose(file_44);

//____________________________________________________________________________
// FASTA  Align from iostream

	Align<String<char>, ArrayGaps> align2;
	fstream strm_1;
	strm_1.open(TEST_PATH "fasta_align.txt", ios_base::in);
	read(strm_1, align2, FastaAlign());
	strm_1.close();
	SEQAN_TASSERT(source(row(align2,0)) == "MKVILLFVLAVFTVFVSSRGIPPEEQSQFLEFQDKFNKKYSHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDEFKNYYLNNKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGPLAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII");
	SEQAN_TASSERT(length(source(row(align2,1))) == 362);
	SEQAN_TASSERT(length(source(row(align2,2))) == 335);

	Align<String<Dna>, ArrayGaps> align22;
	fstream strm_12;
	strm_12.open(TEST_PATH "fasta_align_dna.txt", ios_base::in);
	read(strm_12, align22, FastaAlign());
	strm_12.close();
	SEQAN_TASSERT(source(row(align22,0)) == "CTACGAAAGGTCGTGTCACGATGTCCGCAAGGGATGGCATTGCATAGAGGAATTGATTGCAACCTACGAAA");
	SEQAN_TASSERT(source(row(align22,1)) == "CTTAATGTCCCGCGTACAAGGGATAGCATGTGGCATAGAGGAATAGAATAGCAGCCTACGAAA");

//____________________________________________________________________________
// FASTA Align to C stream

	Align<String<char>, ArrayGaps> align3;
	resize(rows(align3), 2);
	assignSource(row(align3, 0), "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyz");
	insertGaps(row(align3, 0), 2, 3);
	assignSource(row(align3, 1), "xyz");

	String<char> id1("Seq1"); 
	String<char> id2("Seq2"); 
	String<String<char> > str_ids;
	appendValue(str_ids, id1);
	appendValue(str_ids, id2);

	FILE * file_2 = fopen(TEST_PATH "my_fasta_align.txt", "w");
	write(file_2, align3, str_ids, FastaAlign());
	fclose(file_2);

	//test it
	FILE * file_3 = fopen(TEST_PATH "my_fasta_align.txt", "r");
	SEQAN_TASSERT(file_3);

	String<String<char> > str_ids2;
	readIDs(file_3, str_ids2, FastaAlign());
	SEQAN_TASSERT(value(str_ids2,0) == "Seq1");
	SEQAN_TASSERT(value(str_ids2,1) == "Seq2");
	Align<String<char>, ArrayGaps> align4;
	read(file_3, align4, FastaAlign());
	SEQAN_TASSERT(source(row(align4,1)) == "xyz");

	fclose(file_3);

//____________________________________________________________________________
// Dna5 File format output

	Align<String<Dna>, ArrayGaps> align5;
	resize(rows(align5), 2);
	assignSource(row(align5, 0), "aaccggtt");
	assignSource(row(align5, 1), "accgtttt");
	globalAlignment(align5, SimpleScore(), NeedlemanWunsch());


//____________________________________________________________________________
// Test virtual file format object //outdated!!
	/*
	FileFormat<FILE *, String<char>, FastaAlign> fasta_format;
	FileFormat<FILE *, String<char> > & ff = fasta_format;

	//write
	FILE * file_4 = fopen(TEST_PATH "my_fasta_align.txt", "w");

	Align<String<char>, ArrayGaps> align5;
	resize(rows(align5), 2);
	assignSource(row(align5, 0), "abc");
	insertGaps(row(align5, 0), 2, 3);
	assignSource(row(align5, 1), "xyz");

	String<char> id4_1("Seq1"); 
	String<char> id4_2("Seq2"); 
	String<String<char> > str_ids4;
	appendValue(str_ids4, id4_1);
	appendValue(str_ids4, id4_2);

	write(file_4, align5, str_ids4, ff);
	fclose(file_4);
	*/
	// ??? ToDo ???
}


//////////////////////////////////////////////////////////////////////////////

void Test_CGViz() {
//____________________________________________________________________________
// Read FASTA Align from C stream

	FILE * file_1 = fopen(TEST_PATH "fasta_align.txt", "r");
	SEQAN_TASSERT(file_1);
	
	Align<String<char>, ArrayGaps> align;
	read(file_1, align, FastaAlign());
	SEQAN_TASSERT(source(row(align,0)) == "MKVILLFVLAVFTVFVSSRGIPPEEQSQFLEFQDKFNKKYSHEEYLERFEIFKSNLGKIEELNLIAINHKADTKFGVNKFADLSSDEFKNYYLNNKEAIFTDDLPVADYLDDEFINSIPTAFDWRTRGAVTPVKNQGQCGSCWSFSTTGNVEGQHFISQNKLVSLSEQNLVDCDHECMEYEGEEACDEGCNGGLQPNAYNYIIKNGGIQTESSYPYTAETGTQCNFNSANIGAKISNFTMIPKNETVMAGYIVSTGPLAIAADAVEWQFYIGGVFDIPCNPNSLDHGILIVGYSAKNTIFRKNMPYWIVKNSWGADWGEQGYIYLRRGKNTCGVSNFVSTSII");
	SEQAN_TASSERT(length(source(row(align,1))) == 362);
	SEQAN_TASSERT(length(source(row(align,2))) == 335);
	String<String<char> > str_ids1;
	readIDs(file_1, str_ids1, FastaAlign());
	SEQAN_TASSERT(value(str_ids1,0) == "CYS1_DICDI");
	SEQAN_TASSERT(value(str_ids1,1) == "ALEU_HORVU");
	SEQAN_TASSERT(value(str_ids1,2) == "CATH_HUMAN");
	fclose(file_1);


//____________________________________________________________________________
// Output in CGViz format

	FILE * file_2 = fopen(TEST_PATH "my_cgviz_align.txt", "w");
	write(file_2, align, str_ids1, CGViz());
	fclose(file_2);

// To read this file in CGViz:
// java -jar CGViz-1.0beta1.jar -i <data_file>
//
// Reference:
// CGViz - a highly configurable viewer for biomolecular data
// www-ab.informatik.uni-tuebingen.de/software/cgviz/welcome.html
}

//////////////////////////////////////////////////////////////////////////////

void Test_Embl(char * fl_path, char *fl_out_path)
{
	FILE * fl = fopen(fl_path, "rb");
	SEQAN_TASSERT(fl);

	FILE *fl_out = fopen(fl_out_path, "wb");
	SEQAN_TASSERT(fl_out);

	String<char> data;
	String<char> meta;

	readMeta(fl, meta, Embl());
	read(fl, data, Embl());

	SEQAN_TASSERT(length(meta) == 1652)
	SEQAN_TASSERT(length(data) == 281)

	write(fl_out, data, meta, Embl());

	readMeta(fl, meta, Embl());
	read(fl, data, Embl());

	SEQAN_TASSERT(data == "ACGT")

	write(fl_out, data, meta, Embl());

	fclose(fl_out);
	fclose(fl);

	SEQAN_TASSERT(_compareTextFiles(fl_path, fl_out_path))

//____________________________________________________________________________
// test reading data without reading meta before

	fl = fopen(fl_path, "rb");
	SEQAN_TASSERT(fl);

	read(fl, data, Embl());

	SEQAN_TASSERT(length(data) == 281)
	fclose(fl);

}

//////////////////////////////////////////////////////////////////////////////

void Test_Genbank(char * fl_path, char * fl_out_path)
{
	FILE * fl = fopen(fl_path, "rb");
	SEQAN_TASSERT(fl);

	FILE *fl_out = fopen(fl_out_path, "wb");
	SEQAN_TASSERT(fl_out);

	String<char> data;
	String<char> meta;

	readMeta(fl, meta, Genbank());
	read(fl, data, Genbank());

	write(fl_out, data, meta, Genbank());

	SEQAN_TASSERT(length(meta) == 4167)
	SEQAN_TASSERT(length(data) == 5028)

	readMeta(fl, meta, Embl());
	read(fl, data, Genbank());
	write(fl_out, data, meta, Genbank());

	SEQAN_TASSERT(data == "ACGT")

	readMeta(fl, meta, Embl());
	read(fl, data, Genbank());
	write(fl_out, data, meta, Genbank());

	SEQAN_TASSERT(data == "CATAGAT")

	fclose(fl_out);
	fclose(fl);

	SEQAN_TASSERT(_compareTextFiles(fl_path, fl_out_path))
}

//////////////////////////////////////////////////////////////////////////////

//____________________________________________________________________________
//test FileReader iterator
void Test_FileReader_Iterator()
{
	FILE * file_fasta = fopen(TEST_PATH "fasta_crlf.txt", "rb");
	SEQAN_TASSERT(file_fasta);

	Iter<FILE *, FileReader<Fasta> > it_fasta(file_fasta);

	//read the first record "ACGT" from file
	SEQAN_TASSERT(!atEnd(it_fasta))
	SEQAN_TASSERT(value(it_fasta) == 'A');
	goNext(it_fasta);
	SEQAN_TASSERT(!atEnd(it_fasta))
	SEQAN_TASSERT(value(it_fasta) == 'C');
	goNext(it_fasta);
	SEQAN_TASSERT(!atEnd(it_fasta))
	SEQAN_TASSERT(value(it_fasta) == 'G');
	goNext(it_fasta);
	SEQAN_TASSERT(!atEnd(it_fasta))
	SEQAN_TASSERT(value(it_fasta) == 'T');
	goNext(it_fasta);
	SEQAN_TASSERT(atEnd(it_fasta))

	//scan second record = 1065 characters
	goBegin(it_fasta);
	SEQAN_TASSERT(!atEnd(it_fasta))

	unsigned int i = 0;
	while (!atEnd(it_fasta))
	{
		++i;
		goNext(it_fasta);
	}
	SEQAN_TASSERT(i == 1065)


	fclose(file_fasta);
}

//____________________________________________________________________________
//test FileReader string

void Test_FileReader_String()
{
	FILE * file_fasta = fopen(TEST_PATH "fasta_crlf.txt", "rb");
	SEQAN_TASSERT(file_fasta)
	goNext(file_fasta, Fasta());

	String<AminoAcid, FileReader<Fasta> > str(file_fasta);

	SEQAN_TASSERT(length(str) == 1065)
	SEQAN_TASSERT(value(str, 1064) == 'G')

	Iterator<String<AminoAcid, FileReader<Fasta> > >::Type it = begin(str);
	SEQAN_TASSERT(*it == value(str, 0))

	goEnd(it);
	SEQAN_TASSERT(atEnd(it))

	unsigned int pos = length(str);
	do
	{
		--it;
		--pos;
		SEQAN_TASSERT(position(it) == pos)
		SEQAN_TASSERT(*it == value(str, pos))
	} while (!atBegin(it));

	SEQAN_TASSERT(pos == 0);

	while (it < end(str))
	{
		SEQAN_TASSERT(position(it) == pos)
		SEQAN_TASSERT(*it == value(str, pos))
		++pos;
		++it;
	}
	SEQAN_TASSERT(pos == length(str))
	SEQAN_TASSERT(atEnd(it))

	fclose(file_fasta);
}

//____________________________________________________________________________

void Create_Long_Testfile(Fasta)
{
	::std::cout << "create fasta testfile";

	FILE * fl = fopen(TEST_PATH "testfile.txt", "wb");
	fprintf(fl, ">testfile\n"); 
	for (unsigned int i = 0; i < 20; ++i)
	{
		::std::cout << ".";
		for (unsigned int j = 0; j < 500; ++j)
		{
			fprintf(fl, "A.........B.........C.........D.........E.........F.........G.........H.........I.........J.........\n");
		}
		fprintf(fl, "A.........B.........CHERE_IS_THE_PATTERNE.........F.........G.........H.........I.........J.........\n");
		for (unsigned int j = 501; j < 1000; ++j)
		{
			fprintf(fl, "A.........B.........C.........D.........E.........F.........G.........H.........I.........J.........\n");
		}
	}
	fclose(fl);
	::std::cout << "created.\n";
}

void Create_Long_Testfile(Embl)
{
	::std::cout << "create embl testfile";

	FILE * fl = fopen(TEST_PATH "testfile.txt", "wb");
	fprintf(fl, "ID\nSQ\n"); 
	for (unsigned int i = 0; i < 20; ++i)
	{
		::std::cout << ".";
		for (unsigned int j = 0; j < 500; ++j)
		{
			fprintf(fl, "  A.........B.........C.........D.........E.........F.........G.........H.........I.........J.........\n");
		}
		fprintf(fl, "  A.........B.........CHERE_IS_THE_PATTERNE.........F.........G.........H.........I.........J.........\n");
		for (unsigned int j = 501; j < 1000; ++j)
		{
			fprintf(fl, "  A.........B.........C.........D.........E.........F.........G.........H.........I.........J.........\n");
		}
	}
	fclose(fl);
	::std::cout << "created.\n";
}

void Create_Long_Testfile(Genbank)
{
	::std::cout << "create genbank testfile";

	FILE * fl = fopen(TEST_PATH "testfile.txt", "wb");
	fprintf(fl, "LOCUS\nORIGIN\n"); 
	for (unsigned int i = 0; i < 20; ++i)
	{
		::std::cout << ".";
		for (unsigned int j = 0; j < 500; ++j)
		{
			fprintf(fl, "  A.........B.........C.........D.........E.........F.........G.........H.........I.........J.........\n");
		}
		fprintf(fl, "  A.........B.........CHERE_IS_THE_PATTERNE.........F.........G.........H.........I.........J.........\n");
		for (unsigned int j = 501; j < 1000; ++j)
		{
			fprintf(fl, "  A.........B.........C.........D.........E.........F.........G.........H.........I.........J.........\n");
		}
	}
	fclose(fl);
	::std::cout << "created.\n";
}
//____________________________________________________________________________

//test filereader string on long files
template <typename TFormat>
void Test_FileReader_String2()
{
	Create_Long_Testfile(TFormat());

//____________________________________________________________________________
// scan the file using file reader string

	String<char, FileReader<TFormat> > fr(TEST_PATH "testfile.txt");
	Finder<String<char, FileReader<TFormat> > > fnd(fr);
	String<char> ndl = "HERE_IS_THE_PATTERN";
	Pattern<String<char>, Horspool> pat(ndl);

	::std::cout << "find patterns";

	for (unsigned int i = 0; i < 20; ++i)
	{
		bool found = find(fnd, pat);

	::std::cout << ".";
		SEQAN_TASSERT(found);
		SEQAN_TASSERT(!fr.data_scanned);
		SEQAN_TASSERT(position(fnd) == (50021 + i * 100000));
	}

	SEQAN_TASSERT(!find(fnd, pat));
	SEQAN_TASSERT(fr.data_scanned);

	::std::cout << "\n";

//____________________________________________________________________________
// scan the file the old way: read the complete file into memory

	String<char> hay;
	FILE * fl = fopen(TEST_PATH "testfile.txt", "rb");
	read(fl, hay, TFormat());
	fclose(fl);

	Finder<String<char> > fnd2(hay);
	for (unsigned int i = 0; i < 20; ++i)
	{
		bool found = find(fnd2, pat);

	::std::cout << ".";
		SEQAN_TASSERT(found);
		SEQAN_TASSERT(position(fnd2) == (50021 + i * 100000));
	}

}

//____________________________________________________________________________

void Test_FileReader_String3()
{
	FILE * file_fasta = fopen(TEST_PATH "fasta_crlf.txt", "rb");
	SEQAN_TASSERT(file_fasta)

	String<AminoAcid, FileReader<Fasta> > str(file_fasta);
	String<AminoAcid> str2(str);

	SEQAN_TASSERT(str == str2);

	append(str2, str);
	SEQAN_TASSERT(str2 == "ACGTACGT");

	fclose(file_fasta);
}

//____________________________________________________________________________

void Test_FileReader()
{
	Test_FileReader_Iterator();
	Test_FileReader_String();

	Test_FileReader_String3();

	Test_FileReader_String2<Fasta>();
	Test_FileReader_String2<Embl>();
	Test_FileReader_String2<Genbank>();
}

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_Stream();
	Test_CStream();

	Test_Raw();

	Test_Fasta_Read(TEST_PATH "fasta_crlf.txt");
	Test_Fasta_Read(TEST_PATH "fasta_lf.txt");
	Test_Fasta_Read(TEST_PATH "fasta_cr.txt");
	Test_Fasta_Write();

	Test_FastaAlign();
	Test_CGViz();

	Test_Embl(TEST_PATH "embl_crlf.txt", TEST_PATH "test_output.embl.txt");
	Test_Genbank(TEST_PATH "genbank_crlf.txt", TEST_PATH "test_output.genbank.txt");

//	Test_FileReader();

	Test_Embl();
//____________________________________________________________________________

	
	debug::verifyCheckpoints(LIB_PATH "stream.h");
	debug::verifyCheckpoints(LIB_PATH "cstream.h");
	debug::verifyCheckpoints(LIB_PATH "meta.h");
	debug::verifyCheckpoints(LIB_PATH "file_base.h");
	debug::verifyCheckpoints(LIB_PATH "file_format_raw.h");
	debug::verifyCheckpoints(LIB_PATH "file_format_fasta.h");
	debug::verifyCheckpoints(LIB_PATH "file_format_fasta_align.h");
	debug::verifyCheckpoints(LIB_PATH "file_format_embl.h");
	debug::verifyCheckpoints(LIB_PATH "file_format_cgviz.h");
	debug::verifyCheckpoints(LIB_PATH "file_format_guess.h");

	SEQAN_TREPORT("TEST END")

	return 0;
}
