#define SEQAN_DEBUG
#define SEQAN_TEST

#include<seqan/seeds.h>





using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////


void test_global_seed_chain()
{
	typedef Seed<int, SimpleSeed> TSeed;
	String<TSeed> chain;
	String<TSeed> chain2;

	SeedSet<int, SimpleSeed, DefaultScore, void> seedContainer2(1,1);
	addSeed(seedContainer2,1,1,2,3,Single());
	addSeed(seedContainer2,2,5,2,5,Single());
	addSeed(seedContainer2,3,9,2,3,Single());
	addSeed(seedContainer2,5,2,3,4,Single());
	addSeed(seedContainer2,5,8,3,4,Single());
	addSeed(seedContainer2,6,5,1,2,Single());
	addSeed(seedContainer2,8,6,1,2,Single());
	addSeed(seedContainer2,9,1,3,3,Single());
	addSeed(seedContainer2,10,4,2,3,Single());
	addSeed(seedContainer2,10,6,1,2,Single());
	addSeed(seedContainer2,10,8,2,3,Single());

	SEQAN_TASSERT(globalChaining(seedContainer2, chain)==10);
	SEQAN_TASSERT(length(chain) == 4); 

	SEQAN_TASSERT(globalChaining(seedContainer2, chain2, -1 ,13, 13) == -4);

	SEQAN_TASSERT(length(chain2) == 4); 

	String<TSeed> chain3;
	String<TSeed> chain4;

	SeedSet<int, SimpleSeed, DefaultScore, void> seedContainer3(1,1);
	addSeed(seedContainer3,1,1,2,3,Single());

	SEQAN_TASSERT(globalChaining(seedContainer3, chain3)==3);
	SEQAN_TASSERT(length(chain3) == 1); 
}




void Main_GlobalSeedChain(){

	//Fragment<unsigned int, Default()> x();
	
	SEQAN_TREPORT("TEST BEGIN")
	test_global_seed_chain();
	debug::verifyCheckpoints("projects/library/seqan/seeds/global_seed_chain.h");
	SEQAN_TREPORT("TEST END")
}
