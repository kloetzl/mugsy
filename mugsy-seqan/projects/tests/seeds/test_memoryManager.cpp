

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/seeds.h>

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////


void Test_MemoryManager_FreeMemoryPointer()
{
	MemoryManager<int, Block<32>, FreeMemoryPointer > manager;

	SEQAN_TASSERT(capacity(manager) == 0);
	SEQAN_TASSERT(length(manager) == 0);


	Size<MemoryManager<int, Block<32>, FreeMemoryPointer> >::Type id = obtainID(manager);
	
	SEQAN_TASSERT(capacity(manager) == 32);
	SEQAN_TASSERT(length(manager) == 1);

	manager[id] = 4;
	SEQAN_TASSERT(manager[0] == 4);

	obtainID(manager);

	releaseID(manager, 0);
	
	id = obtainID(manager);
	SEQAN_TASSERT(id == 0);


	assignValue(manager, id, 4);
	SEQAN_TASSERT(value(manager,id) == 4);


	clear(manager);
	SEQAN_TASSERT(capacity(manager) == 0);
	

	for(int i = 0; i < 35;i++)
		manager[obtainID(manager)] = i;

	SEQAN_TASSERT(capacity(manager) == 64);
	SEQAN_TASSERT(length(manager) == 35);

	releaseID(manager,6);
	releaseID(manager,11);
	id = obtainID(manager);
	SEQAN_TASSERT(id == 11);
	id =obtainID(manager);

	SEQAN_TASSERT(id == 6);
	
	releaseID(manager,3);
	releaseID(manager,14);
	MemoryManager<int, Block<8>, FreeMemoryPointer> manager2(manager);
	
	id =obtainID(manager2);
	SEQAN_TASSERT(id == 14);
	id =obtainID(manager2);
	SEQAN_TASSERT(id == 3);
	id =obtainID(manager2);
	SEQAN_TASSERT(id == 35);

        typedef MemoryManager<int, Block<8>, FreeMemoryPointer > const TMemoryManager2;
	TMemoryManager2 manager3(manager);
	TMemoryManager2 manager4(manager3);
}



void Test_MemoryManager_FreeMemoryInt()
{
	MemoryManager<int, Block<32>, FreeMemoryInt > manager;

	SEQAN_TASSERT(capacity(manager) == 0);
	SEQAN_TASSERT(length(manager) == 0);


	Size<MemoryManager<int, Block<32>, FreeMemoryInt> >::Type id = obtainID(manager);
	
	SEQAN_TASSERT(capacity(manager) == 32);
	SEQAN_TASSERT(length(manager) == 1);

	manager[id] = 4;
	SEQAN_TASSERT(manager[0] == 4);

	obtainID(manager);

	releaseID(manager, 0);
	
	id = obtainID(manager);
	SEQAN_TASSERT(id == 0);


	assignValue(manager, id, 4);
	SEQAN_TASSERT(value(manager,id) == 4);


	clear(manager);
	SEQAN_TASSERT(capacity(manager) == 0);
	

	for(int i = 0; i < 35;i++)
		manager[obtainID(manager)] = i;
	SEQAN_TASSERT(capacity(manager) == 64);
	SEQAN_TASSERT(length(manager) == 35);

	releaseID(manager,6);
	releaseID(manager,11);
	id = obtainID(manager);
	SEQAN_TASSERT(id == 11);
	id =obtainID(manager);

	SEQAN_TASSERT(id == 6);
	
	releaseID(manager,3);
	releaseID(manager,14);
	MemoryManager<int, Block<8>, FreeMemoryInt> manager2(manager);
	
	id =obtainID(manager2);
	SEQAN_TASSERT(id == 14);
	id =obtainID(manager2);
	SEQAN_TASSERT(id == 3);
	id =obtainID(manager2);
	SEQAN_TASSERT(id == 35);

	MemoryManager<int, Block<8>, FreeMemoryInt > const manager3(manager);
	
	MemoryManager<int, Block<8>, FreeMemoryInt > const manager4(manager3);
        clear(manager);
        clear(manager2);
}

void Main_MemoryManager(){
	SEQAN_TREPORT("TEST BEGIN")
	
	Test_MemoryManager_FreeMemoryPointer();
	Test_MemoryManager_FreeMemoryInt();
	debug::verifyCheckpoints("projects/library/seqan/seeds/memoryManager_base.h");
	debug::verifyCheckpoints("projects/library/seqan/seeds/memoryManager_int.h");
	SEQAN_TREPORT("TEST END")
}
