#ifndef SEQAN_HEADER_TEST_GRAPH_MATCH_REFINEMENT_H
#define SEQAN_HEADER_TEST_GRAPH_MATCH_REFINEMENT_H

#include <seqan/align.h>


using namespace std;
using namespace seqan;

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

int Test_ConvertSequences(String<char> const in_path, String<char> const in_file, String<char> const path, String<char> const file_prefix) {
	typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
	
	// count sequences
	unsigned seqCount = 0;

	ifstream file;
	std::stringstream input;
	input << in_path << in_file;
	file.open(input.str().c_str(), ios_base::in | ios_base::binary);
	if (!file.is_open()) return 0;
	while (!_streamEOF(file)) {
		String<char> id;
		readID(file, id, Fasta());
		std::cout << id << std::endl;
		goNext(file, Fasta());
		++seqCount;
	}
	std::cout << "Number of sequences: " << seqCount << std::endl;

	// import sequences
	file.clear();
	file.seekg(0, ios_base::beg);
	for(unsigned i = 0; (i < seqCount) && !_streamEOF(file); ++i) 
	{
		TString str;
		//String<TraceBack, External<> > trace;
		//open(trace, "D:\\seqan.dat");
		std::stringstream s;
		s << path << file_prefix << i;
		open(str, s.str().c_str());
		read(file, str, Fasta());
	}
    file.close();

	return seqCount;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TVal1, typename TVal2>
inline bool 
Test_ReadSequences(String<char> const path, String<char> const file_prefix, TStringSet& str, TVal1 const start, TVal2 const nseq) {
	for(unsigned i = start; i < start + nseq; ++i) {
		std::stringstream s;
		s << path << file_prefix << i - start;
		bool f = open(str[i], s.str().c_str());
		if (!f) return false;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphMatchRefine() {
	// Sequences
	typedef String<Dna5, External<ExternalConfig<File<>, 64*1024> > > TString;
	typedef StringSet<TString, Owner<> > TStringSet;
	typedef Id<TStringSet>::Type TId;
	typedef Size<TStringSet>::Type TSize;

	// Matches
	typedef String<Fragment<>, External<ExternalConfig<File<>, 64*1024> > > TFragmentString;
	//typedef String<Fragment<>, External<> > TFragmentString;

	// Windows
#ifdef PLATFORM_WINDOWS
	String<char> in_path("Z:\\matches\\");
	String<char> out_path("Z:\\matches\\out\\");
	return;
#else
	// Linux
	String<char> in_path("/home/takifugu2/data/SeqAn/binary/");
	String<char> out_path("/home/takifugu2/data/SeqAn/binary/");
#endif
	

	TSize hSeq = 24;
	TSize wSeq = 24;
	TSize bSeq = 24;
	

	// Convert all sequences only once
	//TSize tmp = 0;
	//tmp = Test_ConvertSequences(in_path, "HUREF6CHROM.fasta",out_path,"H.chr.");
	//if (tmp != hSeq) {
	//	std::cerr << "Did not read all HUREF sequences." << std::endl;
	//	exit(1);
	//}
	//tmp = Test_ConvertSequences(in_path, "WGSACHROM.fasta",out_path,"W.chr.");
	//if (tmp != wSeq) {
	//	std::cerr << "Did not read all WGSA sequences." << std::endl;
	//	exit(1);
	//}
	//tmp = Test_ConvertSequences(in_path, "B36LCCHROM.fasta",out_path,"B.chr.");
	//if (tmp != bSeq) {
	//	std::cerr << "Did not read all B36LC sequences." << std::endl;
	//	exit(1);
	//}

	// Read all the sequences
	TStringSet str;
	resize(str, hSeq + wSeq + bSeq);
	bool f = Test_ReadSequences(out_path,"H.chr.", str, 0, hSeq);
	if (!f) {
		std::cerr << "Error importing HUREF sequences." << std::endl;
		exit(1);
	}
	f = Test_ReadSequences(out_path,"W.chr.", str, hSeq, wSeq);
	if (!f) {
		std::cerr << "Error importing WGSA sequences." << std::endl;
		exit(1);
	}
	f = Test_ReadSequences(out_path,"B.chr.", str, hSeq + wSeq, bSeq);
	if (!f) {
		std::cerr << "Error importing B36LC sequences." << std::endl;
		exit(1);
	}

	// Build a map:
	// SeqId -> Identifier
	typedef std::map<TId, String<char> > TIdToNameMap;
	TIdToNameMap idToName;
	for(TId i = 0; i < length(str); ++i) {
		TSize index = 0;
		std::stringstream s;
		if (i < 24) {
			s << "H";
			index = i;
		}
		else if (i < 48) {
			s << "W";
			index = i - hSeq;
		}
		else {
			s << "B";
			index = i - (hSeq + wSeq);
		}
		s << ":" << index;
		idToName.insert(std::make_pair(i, s.str().c_str()));
	}

	// Just to check that everything worked
	std::cout << "Number of sequences: " << length(str) << std::endl;
	for(TIdToNameMap::const_iterator pos =  idToName.begin(); pos != idToName.end(); ++pos) {
		std::cout << pos->second << ") ";
		for(TSize i=0; ((i<10) && (i <length(str[pos->first])));++i) {
			std::cout << str[pos->first][i];
		}
		std::cout << std::endl;
	}

	// Access the matches
	TFragmentString matches;
	std::stringstream strstream;
	//strstream << out_path << "matchesTest.dat"; // 10 Matches
	//strstream << out_path << "matches1000.dat"; // 2001948 Matches
	strstream << out_path << "matches10000.dat"; // 2111 Matches
	//strstream << out_path << "matches2000.dat"; // 653095 Matches
	//strstream << out_path << "matches500.dat"; // 3999176 
	open(matches, strstream.str().c_str());


	// Convert the matches to an external string
	//for(TSize i = 1; i<4; ++i) {
	//	fstream strm; 
	//	std::stringstream s;
	//	if (i==0 ) s << in_path << "TvsT.atac";
	//	else if (i==1 ) s << in_path << "BvsH.atac";
	//	else if (i==2 ) s << in_path << "BvsW.atac";
	//	else if (i==3 ) s << in_path << "WvsH.atac";
	//	strm.open(s.str().c_str(), ios_base::in);
	//	read(strm, matches, 500, AtacMatches());
	//	strm.close();
	//}

	// Print number of matches
	std::cout << "Number of matches: " << length(matches) << std::endl;
	
	// Re7finement
	typedef Infix<TString>::Type TInfix;
	typedef StringSet<TString, Dependent<> > TAlignmentStringSet;
	typedef Graph<Alignment<TAlignmentStringSet> > TAliGraph;
	typedef VertexDescriptor<TAliGraph>::Type TVD;
	TAlignmentStringSet aliStr;
	for(TSize i = 0; i<length(str); ++i) {
		assignValueById(aliStr, str, i);
	}
	Score<int> score_type = Score<int>(1,-1,-2,0) ;
	TAliGraph ali_graph(aliStr);
	matchRefinement(matches,str,score_type,ali_graph);//,StoreEdges());
	std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
	std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
	//std::cout << ali_graph <<"\n";

	// Print all the matches
	//typedef Iterator<TAliGraph, EdgeIterator>::Type TEdgeIterator;
	//TEdgeIterator it(ali_graph);
	//for(;!atEnd(it);goNext(it)) {
	//	TVD sV = sourceVertex(it);
	//	TVD tV = targetVertex(it);
	//	TId seqId1 = sequenceId(ali_graph,sV);
	//	TId seqId2 = sequenceId(ali_graph,tV);
	//	TSize seqBegin1 = fragmentBegin(ali_graph, sV);
	//	TSize seqBegin2 = fragmentBegin(ali_graph, tV);
	//	TSize len = fragmentLength(ali_graph, sV);
	//	TIdToNameMap::const_iterator pos1 =  idToName.find(seqId1);
	//	TIdToNameMap::const_iterator pos2 =  idToName.find(seqId2);
	//	std::cout << pos1->second << "," << seqBegin1  << "," << len << "," << pos2->second << "," << seqBegin2 << "," << len << std::endl;
	//}

	for(TIdToNameMap::const_iterator pos =  idToName.begin(); pos != idToName.end(); ++pos) {
		close(str[pos->first]);
	}
	close(matches);
}


//produce pairwise alignments (Align object)
template<typename TAlign, typename TSequence, typename TSeqSpec, typename TScore>
void 
getAlignments(String<TAlign> & alis, StringSet<TSequence, TSeqSpec> & seq, TScore & score_type, int & numAlignments, int cutoff)
{

	unsigned int gesamt = 0;

	for(unsigned int i = 0; i < length(seq); ++i)
	{
		for(unsigned int j = i+1; j < length(seq); ++j)
		{
			TAlign ali;
			resize(rows(ali), 2);
			setSource(row(ali, 0), seq[i]);
			setSource(row(ali, 1), seq[j]);

			LocalAlignmentFinder<int> sw_finder = LocalAlignmentFinder<int>(ali);
			
			int score = smithWaterman(ali,sw_finder,score_type,cutoff);
			if(score==0) continue;
			//cout << ali<<"\n";
			//cout <<"Seq "<<i<<" - Seq "<<j<<"\n"<<score<< ali;
			//cout << sourceBeginPosition(row(ali,0)) <<"   ";
			//cout << sourceBeginPosition(row(ali,1)) <<"\n";
			appendValue(alis,ali);
			++gesamt;
			int k = 1;
			while(k<numAlignments)
			{
				score = smithWatermanGetNext(ali,sw_finder,score_type,cutoff);
				if(score==0) break;
				//cout <<score<< ali;
				//cout << sourceBeginPosition(row(ali,0)) <<"   ";
				//cout << sourceBeginPosition(row(ali,1)) <<"\n";
				appendValue(alis,ali);
				++gesamt;
				++k;
			}
		}	
	}

	numAlignments = gesamt;
	resize(alis,numAlignments);


}


void 
Test_RefineAlign(){

	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef Align<TString, ArrayGaps> TAlign;

	int numSequences = 4;

	TStringSet seq_set;


	TString str = "GARFIELDTHELASTFATCAT";
	appendValue(seq_set,str);

	str = "GARFIELDTHEFASTCAT";
	appendValue(seq_set,str);
	
	str = "GARFIELDTHEVERYFASTCAT";
	appendValue(seq_set,str);
	
	str = "THEFATCAT";
	appendValue(seq_set,str);




	int numAlignments = 1;
	int numSequencePairs = 0;
	int cutoff = 3;
	for(int i = 1 ; i < numSequences; ++i) 
		numSequencePairs += i;
	String<TAlign> alis;
	reserve(alis,numSequencePairs*numAlignments);
	Score<int> score_type = Score<int>(1,-1,-2,-2) ;

	getAlignments(alis,seq_set,score_type,numAlignments,cutoff);

	typedef Graph<Alignment<TStringSet> > TAliGraph;
	TAliGraph ali_graph(seq_set);

	//std::cout <<"Number of Segments: "<<length(alis)<<"\n";
	matchRefinement(alis,seq_set,score_type,ali_graph);

	//std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
	//std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
	//std::cout << ali_graph << "\n\n";

	//fstream strmW; // Write the library
	//strmW.open(TEST_PATH "my_testlib1.lib", ios_base::out | ios_base::trunc);
	//write(strmW,ali_graph,TCoffeeLib());
	//strmW.close();

	VertexDescriptor<TAliGraph>::Type vd;

	vd = findVertex(ali_graph,0,0);
	SEQAN_TASSERT(vd == 0)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 8)
	vd = findVertex(ali_graph,0,8);
	SEQAN_TASSERT(vd == 1)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
	vd = findVertex(ali_graph,0,11);
	SEQAN_TASSERT(vd == 2)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
	vd = findVertex(ali_graph,0,14);
	SEQAN_TASSERT(vd == 3)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
	vd = findVertex(ali_graph,0,15);
	SEQAN_TASSERT(vd == 4)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
	vd = findVertex(ali_graph,0,17);
	SEQAN_TASSERT(vd == 5)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
	vd = findVertex(ali_graph,0,18);
	SEQAN_TASSERT(vd == 6)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
	vd = findVertex(ali_graph,0,20);
	SEQAN_TASSERT(vd == 7)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)

	vd = findVertex(ali_graph,1,0);
	SEQAN_TASSERT(vd == 8)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 0)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 8)
	vd = findVertex(ali_graph,1,8);
	SEQAN_TASSERT(vd == 9)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 8)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
	vd = findVertex(ali_graph,1,11);
	SEQAN_TASSERT(vd == 10)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 11)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
	vd = findVertex(ali_graph,1,14);
	SEQAN_TASSERT(vd == 11)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 14)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
	vd = findVertex(ali_graph,1,15);
	SEQAN_TASSERT(vd == 12)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 15)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
	vd = findVertex(ali_graph,1,17);
	SEQAN_TASSERT(vd == 13)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 17)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)

	vd = findVertex(ali_graph,2,0);
	SEQAN_TASSERT(vd == 14)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 0)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 8)
	vd = findVertex(ali_graph,2,8);
	SEQAN_TASSERT(vd == 15)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 8)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
	vd = findVertex(ali_graph,2,11);
	SEQAN_TASSERT(vd == 16)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 11)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 7)
	vd = findVertex(ali_graph,2,18);
	SEQAN_TASSERT(vd == 17)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 18)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
	vd = findVertex(ali_graph,2,19);
	SEQAN_TASSERT(vd == 18)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 19)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
	vd = findVertex(ali_graph,2,21);
	SEQAN_TASSERT(vd == 19)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 21)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)

	vd = findVertex(ali_graph,3,0);
	SEQAN_TASSERT(vd == 20)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 0)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
	vd = findVertex(ali_graph,3,3);
	SEQAN_TASSERT(vd == 21)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 3)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
	vd = findVertex(ali_graph,3,5);
	SEQAN_TASSERT(vd == 22)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 5)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
	vd = findVertex(ali_graph,3,6);
	SEQAN_TASSERT(vd == 23)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 6)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
	vd = findVertex(ali_graph,3,8);
	SEQAN_TASSERT(vd == 24)
	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 8)
	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
	
	SEQAN_TASSERT(findEdge(ali_graph,0,14)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,0,8)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,1,15)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,1,9)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,2,10)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,3,11)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,4,12)!=0)  //doesnt exist if edges with score <= 0 are kicked out
	SEQAN_TASSERT(findEdge(ali_graph,4,21)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,5,22)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,5,13)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,6,23)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,7,24)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,8,14)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,9,20)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,9,15)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,11,22)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,12,23)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,13,24)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,17,22)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,18,23)!=0)
	SEQAN_TASSERT(findEdge(ali_graph,19,24)!=0)


	//clear(ali_graph);
	//assignStringSet(ali_graph,seq_set);

	//matchRefinement(alis,seq_set,score_type,ali_graph,StoreEdges());
	//std::cout << "\nnumEdges: "<<numEdges(ali_graph)<<"\n";
	//std::cout << "\nnumVertices: "<<numVertices(ali_graph)<<"\n";
	//std::cout << ali_graph << "\n\n";


	//fstream strmW1; // Write the library
	//strmW1.open(TEST_PATH "my_testlib2alledges.lib", ios_base::out | ios_base::trunc);
	//write(strmW1,ali_graph,TCoffeeLib());
	//strmW1.close();


}


//produce pairwise alignments (Graph<Alignment>)
//template<typename TAlign, typename TStringSet, typename TScore>
//void 
//getGraphAlignments(String<TAlign> & alis, TStringSet & seq, TScore & score_type, int & numAlignments, int cutoff)
//{
//	typedef StringSet<typename Value<TStringSet>::Type, Dependent<> > TAliStringSet;
//
//	int gesamt = 0;
//
//	for(int i = 0; i < length(seq); ++i)
//	{
//		for(int j = i+1; j < length(seq); ++j)
//		{
//			TAliStringSet str;
//			assignValueById(str, seq[i],positionToId(seq, i));
//			assignValueById(str, seq[j],positionToId(seq, j));
//			TAlign ali_g(str);
//			Value<TScore>::Type score = localAlignment(ali_g, score_type, SmithWatermanClump());
//			if(score==0)
//				continue;
// 			int k = 1;
//			while(k<numAlignments)
//			{
//				score = localAlignment(ali_g, score_type, SmithWatermanClump());
//				if(score==0) k = numAlignments;
//				else ++k;
//			}
//			appendValue(alis,ali_g);
//			cout << ali_g <<"\n";
//			++gesamt;
//		}	
//	}
//
//	numAlignments = gesamt;
//	resize(alis,numAlignments);
//
//
//}
//
//
//
//
//void 
//Test_RefineAlignGraph(){
//
//	typedef String<char> TString;
//	typedef StringSet<TString> TStringSet;
//	//typedef Align<typename Reference<TStringSet>::Type, ArrayGaps> TAlign;
//	//typedef Graph<Alignment<TStringSet, unsigned int> > TAlign;
//	typedef Graph<Alignment<StringSet<TString,Dependent<> >, unsigned int> > TAlign;
//
//
//	TStringSet seq_set;
//
//
//	TString str = "GARFIELDTHELASTFATCAT";
//	//appendValue(seq_set,str);
//	assignValueById(seq_set,str);
//
//	str = "GARFIELDTHEFASTCAT";
//	//appendValue(seq_set,str);
//	assignValueById(seq_set,str);
//
//	str = "GARFIELDTHEVERYFASTCAT";
//	//appendValue(seq_set,str);
//	assignValueById(seq_set,str);
//	
//	str = "THEFATCAT";
//	//appendValue(seq_set,str);
//	assignValueById(seq_set,str);
//
//	int numSequences = length(seq_set);
//
//
//	int numAlignments = 2;
//	int numSequencePairs = 0;
//	int cutoff = 3;
//	for(int i = 1 ; i < numSequences; ++i) 
//		numSequencePairs += i;
//	String<TAlign> alis;
//	reserve(alis,numSequencePairs*numAlignments);
//	Score<int> score_type = Score<int>(1,-1,-2,-2) ;
//
//	getGraphAlignments(alis,seq_set,score_type,numAlignments,cutoff);
//
//	typedef Graph<Alignment<TStringSet> > TAliGraph;
//	//TAliGraph ali_graph;
//	TAliGraph ali_graph(seq_set);
//
//	matchRefinement(alis,seq_set,score_type,ali_graph);
//
//	cout << ali_graph << "\n";
//	VertexDescriptor<TAliGraph>::Type vd;
//
//	vd = findVertex(ali_graph,0,0);
//	SEQAN_TASSERT(vd == 0)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 8)
//	vd = findVertex(ali_graph,0,8);
//	SEQAN_TASSERT(vd == 1)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
//	vd = findVertex(ali_graph,0,11);
//	SEQAN_TASSERT(vd == 2)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
//	vd = findVertex(ali_graph,0,14);
//	SEQAN_TASSERT(vd == 3)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//	vd = findVertex(ali_graph,0,15);
//	SEQAN_TASSERT(vd == 4)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
//	vd = findVertex(ali_graph,0,17);
//	SEQAN_TASSERT(vd == 5)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//	vd = findVertex(ali_graph,0,18);
//	SEQAN_TASSERT(vd == 6)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
//	vd = findVertex(ali_graph,0,20);
//	SEQAN_TASSERT(vd == 7)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//
//	vd = findVertex(ali_graph,1,0);
//	SEQAN_TASSERT(vd == 8)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 0)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 8)
//	vd = findVertex(ali_graph,1,8);
//	SEQAN_TASSERT(vd == 9)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 8)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
//	vd = findVertex(ali_graph,1,11);
//	SEQAN_TASSERT(vd == 10)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 11)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
//	vd = findVertex(ali_graph,1,14);
//	SEQAN_TASSERT(vd == 11)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 14)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//	vd = findVertex(ali_graph,1,15);
//	SEQAN_TASSERT(vd == 12)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 15)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
//	vd = findVertex(ali_graph,1,17);
//	SEQAN_TASSERT(vd == 13)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 17)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//
//	vd = findVertex(ali_graph,2,0);
//	SEQAN_TASSERT(vd == 14)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 0)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 8)
//	vd = findVertex(ali_graph,2,8);
//	SEQAN_TASSERT(vd == 15)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 8)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
//	vd = findVertex(ali_graph,2,11);
//	SEQAN_TASSERT(vd == 16)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 11)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 7)
//	vd = findVertex(ali_graph,2,18);
//	SEQAN_TASSERT(vd == 17)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 18)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//	vd = findVertex(ali_graph,2,19);
//	SEQAN_TASSERT(vd == 18)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 19)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
//	vd = findVertex(ali_graph,2,21);
//	SEQAN_TASSERT(vd == 19)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 21)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//
//	vd = findVertex(ali_graph,3,0);
//	SEQAN_TASSERT(vd == 20)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 0)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 3)
//	vd = findVertex(ali_graph,3,3);
//	SEQAN_TASSERT(vd == 21)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 3)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
//	vd = findVertex(ali_graph,3,5);
//	SEQAN_TASSERT(vd == 22)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 5)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//	vd = findVertex(ali_graph,3,6);
//	SEQAN_TASSERT(vd == 23)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 6)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 2)
//	vd = findVertex(ali_graph,3,8);
//	SEQAN_TASSERT(vd == 24)
//	SEQAN_TASSERT(fragmentBegin(ali_graph,vd) == 8)
//	SEQAN_TASSERT(fragmentLength(ali_graph,vd) == 1)
//	
//	SEQAN_TASSERT(findEdge(ali_graph,0,14)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,0,8)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,1,15)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,1,9)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,2,10)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,3,11)!=0)
////	SEQAN_TASSERT(findEdge(ali_graph,4,12)!=0)  //doesnt exist if edges with score <= 0 are kicked out
//	SEQAN_TASSERT(findEdge(ali_graph,4,21)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,5,22)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,5,13)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,6,23)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,7,24)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,8,14)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,9,20)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,9,15)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,11,22)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,12,23)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,13,24)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,17,22)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,18,23)!=0)
//	SEQAN_TASSERT(findEdge(ali_graph,19,24)!=0)
//
//
//}


void Test_Problem() 
{
	typedef String<char> TString;
	typedef StringSet<TString> TStringSet;
	//typedef StringSet<TString, Dependent<> > TAlignmentStringSet;
	typedef Graph<Alignment<TStringSet> > TAlign;
	typedef Fragment<> TFragment;
	typedef String<TFragment> TFragString;


	TString str1 = "RKNLVQFGVGEKNGSVRWVMNALGVKDDWLLVPSHAYKFEKDYEMMEFYFNRGGTYYSISAGNVVIQSLDVGFQDVVLMKVPTIPKFRDITQHFIKKGDVPRALNRLATLVTTVNGTPMLISEGPLKMEEKATYVHKKNDGTTVDLTVDQAWRGKGEGLPGMCGGALVSSNQSIQNAILGIHVAGGNSILVAKLVTQE";
	TString str2 = "";
	TString str3 = "";
	TString str4 = "IAGGEAITTGGSRCSLGFNVVAHALTAGHCTNISAWSIGTRTGTSFNNDYGIIRHSNPAAADGRVYLYQDITTAGNAFVGQAVQRSGSTTGLRSGSVTGLNATVNYGSSGIVYGMIQTNVCAGDSGGSLFAGSTALGLTSGGSGNCRTGGTTFYQPVT";
	TString str5 = "";

	TStringSet strSet;
	assignValueById(strSet,str1,0);
	assignValueById(strSet,str2,1);
	assignValueById(strSet,str3,2);
	assignValueById(strSet,str4,3);
	assignValueById(strSet,str5,4);
	//cout << length(strSet)<<"\n";
	//cout << idToPosition(strSet,0) <<"\n";
	//cout << idToPosition(strSet,3) <<"\n";
	//cout << strSet[0] <<"\n";
	//cout << strSet[1] <<"\n";
	//cout << strSet[2] <<"\n";
	//cout << strSet[3] <<"\n";
	//cout << positionToId(strSet,0) <<"\n";
	//cout << positionToId(strSet,1) <<"\n";


	TFragment frag(0,28,3,35,18);

//	typedef String<Fragment<>, External<ExternalConfig<File<>, 64*1024> > > TFragString;

	TFragString alis;
	//resize(alis,1);
	//alis[0] = frag;
	//appendValue(alis, (TFragment) frag);

	//TAlign outGraph(strSet);
	//matchRefinement(alis,strSet,outGraph);

}


void Test_GraphMatchRefinement() 
{

//	Test_Problem();

	//test for refinement on Align<TSource,TSpec>
//	Test_RefineAlign();
	
	//test for refinement on Graph<Alignment<> >
//	Test_RefineAlignGraph();

	//test for refinement on Fragment<>
//	Test_GraphMatchRefine();
}

}

#endif

