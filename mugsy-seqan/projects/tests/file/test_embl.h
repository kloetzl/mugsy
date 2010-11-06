#ifndef SEQAN_HEADER_TEST_EMBL_H
#define SEQAN_HEADER_TEST_EMBL_H

#define TEST_PATH "projects/tests/file/"
#define LIB_PATH "projects/library/seqan/file/"

using namespace std;

namespace SEQAN_NAMESPACE_MAIN
{

	

	
//////////////////////////////////////////////////////////////////////////////
void Test_EmblOnFile() {


	std::fstream strm; 
	strm.open(TEST_PATH "takifugu_scl_embl.txt", ios_base::in | ios_base::binary);

	String<char> line;
	String<char> feature_line;

	readLineType(strm, feature_line, "FT", Embl());
	//cout << feature_line << "\n";

	int count = 0;
	int next_pos = readFeature(feature_line, 0, line, "exon", Embl());
	while(next_pos != 0)
	{
		++count;
	//  cout << line << "\n";
		next_pos = readFeature(feature_line, next_pos, line, "exon", Embl());
	}

	SEQAN_TASSERT(count == 3)


}


void Test_EmblOnMeta() {


	std::fstream strm; 
	strm.open(TEST_PATH "takifugu_scl_embl.txt", ios_base::in | ios_base::binary);

	String<char> line;

	String<char> feature_line;
	String<char> meta;
	readMeta(strm,meta,Embl());
	
	readLineType(meta, line, "KW", Embl());
	SEQAN_TASSERT(line == "SCL gene.")

	readLineType(meta, line, "RX", Embl());
	SEQAN_TASSERT(infix(line,0,28) == "DOI; 10.1073/pnas.101532998.")
	SEQAN_TASSERT(length(line)==46 ||length(line)==47)

	clear(line);
	readLineType(meta, feature_line, "FT", Embl());
	//cout << feature_line << "\n";



	int count = 0;
	int next_pos = readFeature(feature_line, 0, line, "CDS", Embl());
	while(next_pos != 0)
	{
		++count;
	//  cout << line << "\n";
		next_pos = readFeature(feature_line, next_pos, line, "CDS", Embl());
	}

	SEQAN_TASSERT(count == 1)


}


void Test_Embl() 
{
	Test_EmblOnFile();
	Test_EmblOnMeta();
}


}

#endif

