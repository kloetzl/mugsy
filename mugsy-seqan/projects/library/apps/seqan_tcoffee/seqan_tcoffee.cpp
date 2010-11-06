/*==========================================================================
               SeqAn - The Library for Sequence Analysis
                         http://www.seqan.de 
============================================================================
Copyright (C) 2007

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
==========================================================================*/


#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include "rna_alphabet.h"
#include <seqan/modifier.h>
#include <seqan/misc/misc_cmdparser.h>

#include <iostream>

using namespace seqan;


//////////////////////////////////////////////////////////////////////////////////

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 4637 $";
	addVersionLine(parser, "Version 1.11 (30. July 2009) Revision: " + rev.substr(11, 4) + "");
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TSeqSet, typename TNameSet>
bool _loadSequences(TSeqSet& sequences, 
					TNameSet& fastaIDs,
					const char *fileName)
{
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);
	unsigned seqCount = length(multiFasta);
	resize(sequences, seqCount, Exact());
	resize(fastaIDs, seqCount, Exact());
	for(unsigned i = 0; i < seqCount; ++i) 
	{
		assignSeqId(fastaIDs[i], multiFasta[i], format);
		assignSeq(sequences[i], multiFasta[i], format);
	}
	return (seqCount > 0);
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
customizedMsaAlignment(MsaOptions<TAlphabet, TScore> const& msaOpt) {
	typedef String<TAlphabet> TSequence;
	StringSet<TSequence, Owner<> > sequenceSet;
	StringSet<String<char> > sequenceNames;
	_loadSequences(sequenceSet, sequenceNames, msaOpt.seqfile.c_str());
	
	// Alignment of the sequences
	Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > gAlign;
	
	// MSA
	globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
		
	// Alignment output
	if (msaOpt.outputFormat == 0) {
		FILE* strmWrite = fopen(msaOpt.outfile.c_str(), "w");
		write(strmWrite, gAlign, sequenceNames, FastaFormat());
		fclose(strmWrite);
	} else if (msaOpt.outputFormat == 1) {
		FILE* strmWrite = fopen(msaOpt.outfile.c_str(), "w");
		write(strmWrite, gAlign, sequenceNames, MsfFormat());
		fclose(strmWrite);
	}

}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, TScore>&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, TScore>&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc msc) {
	msaOpt.sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc mmsc) {
	msaOpt.sc.data_mismatch = mmsc;
}


//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
_initMsaParams(CommandLineParser& parser, TScore& scMat) {
	
	// Msa configuration
	MsaOptions<TAlphabet, TScore> msaOpt;
	
	// Set main options
	getOptionValueLong(parser, "seq", msaOpt.seqfile);
	getOptionValueLong(parser, "outfile", msaOpt.outfile);
	String<char> optionVal;
	getOptionValueLong(parser, "format", optionVal);
	if (optionVal == "fasta") msaOpt.outputFormat = 0;
	else if (optionVal == "msf") msaOpt.outputFormat = 1;

	// Set segment match generation options
	::std::string tmpVal;
	getOptionValueLong(parser, "method", tmpVal);
	unsigned int beg = 0;
	for(unsigned int i = 0; i<tmpVal.length(); ++i) {
		if (tmpVal[i] == ',') {
			if (tmpVal.substr(beg, i - beg) == "global") appendValue(msaOpt.method, 0);
			else if (tmpVal.substr(beg, i - beg) == "local") appendValue(msaOpt.method, 1);
			else if (tmpVal.substr(beg, i - beg) == "overlap") appendValue(msaOpt.method, 2);
			else if (tmpVal.substr(beg, i - beg) == "lcs") appendValue(msaOpt.method, 3);
			beg = i + 1;
		}
	}
	if (beg != tmpVal.length()) {
			if (tmpVal.substr(beg, tmpVal.length() - beg) == "global") appendValue(msaOpt.method, 0);
			else if (tmpVal.substr(beg, tmpVal.length() - beg) == "local") appendValue(msaOpt.method, 1);
			else if (tmpVal.substr(beg, tmpVal.length() - beg) == "overlap") appendValue(msaOpt.method, 2);
			else if (tmpVal.substr(beg, tmpVal.length() - beg) == "lcs") appendValue(msaOpt.method, 3);
	}
	getOptionValueLong(parser, "blast", tmpVal);
	beg = 0;
	for(unsigned int i = 0; i<tmpVal.length(); ++i) {
		if (tmpVal[i] == ',') {
			appendValue(msaOpt.blastfiles, tmpVal.substr(beg, i - beg));
			beg = i + 1;
		}
	}
	if (beg != tmpVal.length())
		appendValue(msaOpt.blastfiles, tmpVal.substr(beg, tmpVal.length() - beg));
	getOptionValueLong(parser, "mummer", tmpVal);
	beg = 0;
	for(unsigned int i = 0; i<tmpVal.length(); ++i) {
		if (tmpVal[i] == ',') {
			appendValue(msaOpt.mummerfiles, tmpVal.substr(beg, i - beg));
			beg = i + 1;
		}
	}
	if (beg != tmpVal.length())
		appendValue(msaOpt.mummerfiles, tmpVal.substr(beg, tmpVal.length() - beg));	
	getOptionValueLong(parser, "lib", tmpVal);
	beg = 0;
	for(unsigned int i = 0; i<tmpVal.length(); ++i) {
		if (tmpVal[i] == ',') {
			appendValue(msaOpt.libfiles, tmpVal.substr(beg, i - beg));
			beg = i + 1;
		}
	}
	if (beg != tmpVal.length())
		appendValue(msaOpt.libfiles, tmpVal.substr(beg, tmpVal.length() - beg));	
	getOptionValueLong(parser, "aln", tmpVal);
	beg = 0;
	for(unsigned int i = 0; i<tmpVal.length(); ++i) {
		if (tmpVal[i] == ',') {
			appendValue(msaOpt.alnfiles, tmpVal.substr(beg, i - beg));
			beg = i + 1;
		}
	}
	if (beg != tmpVal.length())
		appendValue(msaOpt.alnfiles, tmpVal.substr(beg, tmpVal.length() - beg));

	// Set scoring options
	msaOpt.sc = scMat;
	getOptionValueLong(parser, "gop", msaOpt.sc.data_gap_open);
	getOptionValueLong(parser, "gex", msaOpt.sc.data_gap_extend);
	int msc = 0;
	getOptionValueLong(parser, "msc", msc);
	_setMatchScore(msaOpt, msc);
	int mmsc = 0;
	getOptionValueLong(parser, "mmsc", mmsc);
	_setMismatchScore(msaOpt, mmsc);

	// Set guide tree options
	getOptionValueLong(parser, "usetree", msaOpt.treefile);
	getOptionValueLong(parser, "build", optionVal);
	if (optionVal == "nj") msaOpt.build = 0;
	else if (optionVal == "min") msaOpt.build = 1;
	else if (optionVal == "max") msaOpt.build = 2;
	else if (optionVal == "avg") msaOpt.build = 3;
	else if (optionVal == "wavg") msaOpt.build = 4;

	// Set alignment evaluation	options
	getOptionValueLong(parser, "infile", msaOpt.infile);

	// Check if any segment-match generation procedure is selected, otherwise set the default
	if ((empty(msaOpt.blastfiles)) && (empty(msaOpt.mummerfiles)) && (empty(msaOpt.libfiles)) && (empty(msaOpt.alnfiles)) && (empty(msaOpt.method))) {
		appendValue(msaOpt.method, 0);
		appendValue(msaOpt.method, 1);
	}

	// Evaluation mode?
	if (isSetLong(parser, "infile")) {
		evaluateAlignment(msaOpt);
	} else { // or alignment mode?
		if (!isSetLong(parser, "seq")) { 
			shortHelp(parser, std::cerr);	// print short help and exit
			exit(0);
		}
		customizedMsaAlignment(msaOpt);
	}
}


//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(CommandLineParser& parser, Dna5 const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initMsaParams<Dna5>(parser, sc);
	} else {
		Score<int> sc;
		_initMsaParams<Dna5>(parser, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(CommandLineParser& parser, Rna5 const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initMsaParams<Rna5>(parser, sc);
	} else {
		Score<int> sc;
		_initMsaParams<Rna5>(parser, sc);
	}
}

//////////////////////////////////////////////////////////////////////////////////

inline void
_initScoreMatrix(CommandLineParser& parser, AminoAcid const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initMsaParams<AminoAcid>(parser, sc);
	} else {
		Blosum62 sc;
		_initMsaParams<AminoAcid>(parser, sc);
	}
}



//////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char *argv[]) {
	// Command line parsing
	CommandLineParser parser;
	_addVersion(parser);
	
	addTitleLine(parser, "*************************************************");
	addTitleLine(parser, "* Multiple sequence alignment - SeqAn::T-Coffee *");
	addTitleLine(parser, "* (c) Copyright 2009 by Tobias Rausch           *");
	addTitleLine(parser, "*************************************************");

	addUsageLine(parser, "-s <FASTA sequence file> [Options]");

	addSection(parser, "Main Options:");
	addOption(parser, addArgumentText(CommandLineOption("s", "seq", "file with sequences", OptionType::String), "<FASTA Sequence File>"));
	addOption(parser, addArgumentText(CommandLineOption("a", "alphabet", "sequence alphabet", (int)OptionType::String, "protein"), "[protein | dna | rna]"));
	addOption(parser, addArgumentText(CommandLineOption("o", "outfile", "output filename", (int)OptionType::String, "out.fasta"), "<Filename>"));
	addOption(parser, addArgumentText(CommandLineOption("f", "format", "output format", (int)OptionType::String, "fasta"), "[fasta | msf]"));
	
	addSection(parser, "Segment Match Generation Options:");
	addOption(parser, CommandLineOption("m", "method", "list of match generation methods", OptionType::String));
	addHelpLine(parser, "global = Global alignments");
	addHelpLine(parser, "local = Local alignments");
	addHelpLine(parser, "overlap = Overlap alignments");
	addHelpLine(parser, "lcs = Longest common subsequence");
	addHelpLine(parser, "Default: global,local");
	addHelpLine(parser, "/*No spaces in-between.*/");
	addOption(parser, addArgumentText(CommandLineOption("bl", "blast", "list of BLAST match files", OptionType::String), "<File1>,<File2>,..."));
	addOption(parser, addArgumentText(CommandLineOption("mu", "mummer", "list of MUMmer match files", OptionType::String), "<File1>,<File2>,..."));
	addOption(parser, addArgumentText(CommandLineOption("al", "aln", "list of FASTA align files", OptionType::String), "<File1>,<File2>,..."));
	addOption(parser, addArgumentText(CommandLineOption("li", "lib", "list of T-Coffee libraries", OptionType::String), "<File1>,<File2>,..."));

	addSection(parser, "Scoring Options:");
	addOption(parser, addArgumentText(CommandLineOption("g", "gop", "gap open penalty", (int)OptionType::Int, -13), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("e", "gex", "gap extension penalty", (int)OptionType::Int, -1), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("ma", "matrix", "score matrix", (int)OptionType::String, "Blosum62"), "<Matrix file>"));
	addOption(parser, addArgumentText(CommandLineOption("ms", "msc", "match score", (int)OptionType::Int, 5), "<Int>"));
	addOption(parser, addArgumentText(CommandLineOption("mm", "mmsc", "mismatch penalty", (int)OptionType::Int, -4), "<Int>"));
	
	addSection(parser, "Guide Tree Options:");
	addOption(parser, addArgumentText(CommandLineOption("u", "usetree", "tree filename", OptionType::String), "<Newick guide tree>"));
	addOption(parser, addArgumentText(CommandLineOption("b", "build", "tree building method", (int)OptionType::String, "nj"), "[nj, min, max, avg, wavg]"));
	addHelpLine(parser, "nj = Neighbor-joining");
	addHelpLine(parser, "min = UPGMA single linkage");
	addHelpLine(parser, "max = UPGMA complete linkage");
	addHelpLine(parser, "avg = UPGMA average linkage");
	addHelpLine(parser, "wavg = UPGMA weighted average linkage");
	addHelpLine(parser, "/*Neighbor-joining creates an");
	addHelpLine(parser, "  unrooted tree. We root that tree");
	addHelpLine(parser, "  at the last joined pair.*/");

	// Alignment evaluation	
	addSection(parser, "Alignment Evaluation Options:");
	addOption(parser, addArgumentText(CommandLineOption("i", "infile", "alignment file", OptionType::String), "<FASTA alignment file>"));

	if (argc == 1)
	{
		shortHelp(parser, std::cerr);	// print short help and exit
		return 0;
	}

	if (!parse(parser, argc, argv, ::std::cerr)) return 1;
	if (isSetLong(parser, "help") || isSetLong(parser, "version")) return 0;	// print help or version and exit


	// Basic command line options
	String<char> alphabet;
	getOptionValueLong(parser, "alphabet", alphabet);
	
	// Initialize scoring matrices
	if (alphabet == "dna") _initScoreMatrix(parser, Dna5());
	else if (alphabet == "rna") _initScoreMatrix(parser, Rna5());
	else _initScoreMatrix(parser, AminoAcid());

	return 0;
}
