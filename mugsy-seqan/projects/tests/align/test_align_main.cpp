#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/basic.h>

void Main_TestGaps(); //test_align_gaps.cpp
void Main_TestAlign(); //test_align_align.cpp
void Main_TestLocalAlign();//test_align_local.cpp
void Main_TestMyers(); // test_align_myers.cpp

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Main_TestLocalAlign();
	Main_TestGaps();
	Main_TestAlign();
	Main_TestMyers();

	SEQAN_TREPORT("TEST END");

	return 0;
}

