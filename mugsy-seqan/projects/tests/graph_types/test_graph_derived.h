#ifndef SEQAN_HEADER_TEST_GRAPH_DERIVED_H
#define SEQAN_HEADER_TEST_GRAPH_DERIVED_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

void Test_Oracle() {
	Graph<Automaton<char> > g;
	createOracleOnReverse(g,"announce");
	SEQAN_TASSERT(parseString(g, 0, "e") == 1)
	SEQAN_TASSERT(parseString(g, 0, "ec") == 2)
	SEQAN_TASSERT(parseString(g, 0, "n") == 3)
	SEQAN_TASSERT(parseString(g, 0, "a") == 8)
	SEQAN_TASSERT(parseString(g, 0, "nn") == 7)

	Graph<Automaton<Dna> > g2;
	createOracle(g2,"ATATA");
	SEQAN_TASSERT(parseString(g2, 0, "A") == 1)
	SEQAN_TASSERT(parseString(g2, 0, "T") == 2)
	SEQAN_TASSERT(parseString(g2, 0, "AT") == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Trie() {
	Graph<Automaton<char> > g;
	String<String<unsigned int> > pos;
	String<String<char> > keywords;
	appendValue(keywords, String<char>("announce"));
	appendValue(keywords, String<char>("annual"));
	appendValue(keywords, String<char>("annually"));
	createTrie(g,pos,keywords);

	SEQAN_TASSERT(parseString(g, 0, "a") == 1)
	SEQAN_TASSERT(parseString(g, 0, "an") == 2)
	SEQAN_TASSERT(parseString(g, 0, "ann") == 3)
	SEQAN_TASSERT(parseString(g, 0, "anno") == 4)
	SEQAN_TASSERT(parseString(g, 0, "annu") == 9)
	SEQAN_TASSERT(getProperty(pos, 11) == (unsigned int) 1) // In vertex 11 keyword 1 ends
	SEQAN_TASSERT(getProperty(pos, 13) == (unsigned int) 2)
	SEQAN_TASSERT(getProperty(pos, 8) == (unsigned int) 0)

	clear(g);
	clear(pos);
	createTrieOnReverse(g,pos,keywords);

	SEQAN_TASSERT(parseString(g, 0, "e") == 1)
	SEQAN_TASSERT(parseString(g, 0, "l") == 9)
	SEQAN_TASSERT(parseString(g, 0, "y") == 15)
	SEQAN_TASSERT(parseString(g, 0, "ec") == 2)
	SEQAN_TASSERT(getProperty(pos, 8) == (unsigned int) 0) // In vertex 8 keyword 0 ends
	SEQAN_TASSERT(getProperty(pos, 14) == (unsigned int) 1)
	SEQAN_TASSERT(getProperty(pos, 22) == (unsigned int) 2)

	Graph<Automaton<Dna> > gDna;
	clear(pos);
	String<String<Dna> > keyw;
	appendValue(keyw, String<Dna>("ATATATA"));
	appendValue(keyw, String<Dna>("TATAT"));
	appendValue(keyw, String<Dna>("ACGATAT"));
	createTrie(gDna,pos,keyw);

	SEQAN_TASSERT(parseString(gDna, 0, "A") == 1)
	SEQAN_TASSERT(parseString(gDna, 0, "T") == 8)
	SEQAN_TASSERT(parseString(gDna, 0, "AT") == 2)
	SEQAN_TASSERT(parseString(gDna, 0, "AC") == 13)
	SEQAN_TASSERT(getProperty(pos, 7) == (unsigned int) 0) // In vertex 7 keyword 0 ends
	SEQAN_TASSERT(getProperty(pos, 18) == (unsigned int) 2)
	SEQAN_TASSERT(getProperty(pos, 12) == (unsigned int) 1)


	//createSuffixTrie
	clear(g);
	clear(pos);
	char * str = (char*) "ABABBA";
	char * strend = end(str);
	char * it;
	typedef VertexDescriptor<Graph<Automaton<char> > >::Type TVertexDescriptor;
	TVertexDescriptor v;

	createSuffixTrie(g, pos, str);
	for (unsigned int i = 0; i < length(str); ++i)
	{
		v = parseString(g, 0, it=str+i, strend);
		SEQAN_TASSERT(it == strend)
		SEQAN_TASSERT(getProperty(pos, v) == i)
	}

	SEQAN_TASSERT(canParseString(g, "ABBA"))
	SEQAN_TASSERT(canParseString(g, "BAB"))

	SEQAN_TASSERT(!canParseString(g, "AA"))
	SEQAN_TASSERT(!canParseString(g, "BAC"))
	SEQAN_TASSERT(!canParseString(g, "C"))
	SEQAN_TASSERT(canParseString(g, ""))

}

//////////////////////////////////////////////////////////////////////////////

void Test_SetOracle() 
{
	Graph<Automaton<char> > g;
	String<String<unsigned int> > pos;
	String<String<char> > keywords;
	appendValue(keywords, String<char>("announce"));
	appendValue(keywords, String<char>("annual"));
	appendValue(keywords, String<char>("annually"));
	createSetOracle(g,pos,keywords);

	for (unsigned int i = 0; i < length(keywords); ++i)
	{
		String<char> & str = keywords[i];
		for (unsigned int j = 0; j < length(str); ++j)
		{
			SEQAN_TASSERT(canParseString(g, prefix(str, i)));
		}
	}

	SEQAN_TASSERT(!canParseString(g, "d"))
	SEQAN_TASSERT(!canParseString(g, "annly"))
	SEQAN_TASSERT(canParseString(g, ""))

}
//////////////////////////////////////////////////////////////////////////////

void Test_GraphDerivedTypes() {
	Test_Oracle();
	Test_Trie();
	Test_SetOracle();

	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_oracle.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_trie.h");
}

}

#endif

