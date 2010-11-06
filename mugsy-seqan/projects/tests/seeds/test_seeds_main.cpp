#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_TEST

#include <seqan/seeds.h>

void Main_BandedAlign();
void Main_GlobalSeedChain();
void Main_MemoryManager();
void Main_Seeds();
void Main_SeedSet();

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Main_BandedAlign();
	Main_GlobalSeedChain();
	Main_MemoryManager();
	Main_Seeds();
	Main_SeedSet();


	SEQAN_TREPORT("TEST END");

	return 0;
}
