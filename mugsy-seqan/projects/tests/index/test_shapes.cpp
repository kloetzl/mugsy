#include <iostream>
#include <seqan/index.h>

using namespace seqan;
using namespace std;

template <typename TShape1, typename TShape2>
bool testShape(TShape1 shape1, TShape2 shape2, bool dump)
{
	bool result = true;

	if (dump) cout << endl;
	              // 012345678901234 len=15
	DnaString dna = "CGGTACGTAAGTTAG";
	DnaString dna1, dna2;

	TShape1 shape1b(shape1);
	TShape2 shape2b(shape2);

	result &= (length(shape1) == length(shape2));
	result &= (weight(shape1) == weight(shape2));
	
	Iterator<DnaString>::Type it = begin(dna);	
	unsigned H1b, H2b;
	for (int i = length(dna); i >= 0; --i)
	{
		{
			unsigned H1 = hash(shape1, it, i);
			unsigned H2 = hash(shape2, it, i);
			if (i >= (int)length(shape1) && i != (int)length(dna))
			{
				H1b = hashNext(shape1b, it);
				H2b = hashNext(shape2b, it);
			} else {
				H1b = hash(shape1b, it, i);
				H2b = hash(shape2b, it, i);
			}

			unhash(dna1, H1, weight(shape1));
			unhash(dna2, H2, weight(shape2));
			if (dump) cout << dec << i << "\t" << hex << H1 << " " << H2 << ' ' << H1b << ' ' << H2b <<"\t" << dna1 << " " << dna2 << "\t" << (H1 == H2 && H1 == H1b && H1b == H2b);

			result &= (H1 == H2 && H1 == H1b && H1b == H2b);

			if (i >= (int)length(shape1)) 
			{
				unsigned H1c = hash(shape1, it);
				unsigned H2c = hash(shape2, it);
				if (dump) cout << " " << (H1c == H2c && H1 == H1c);
			}
		}
		if (dump) cout << "\t";
		{
			unsigned H1 = hashUpper(shape1, it, i);
			unsigned H2 = hashUpper(shape2, it, i);
			
			unhash(dna1, H1, weight(shape1));
			unhash(dna2, H2, weight(shape2));
			if (dump) cout << dec << i << "\t" << hex << H1 << " " << H2 << "\t" << dna1 << " " << dna2 << "\t" << (H1 == H2);

			result &= (H1 == H2);
		}
		if (dump) cout << endl;
		++it;
	}
	return result;
}

void testShapes() 
{
	Shape<Dna, SimpleShape > shapeA(6);
	SEQAN_TASSERT(testShape(shapeA, Shape<Dna, UngappedShape<6> >(), true));

	                   // 012345678  len=9
	CharString pattern = "11100110100";
	Shape<Dna, GenericShape> shapeB(pattern);
	SEQAN_TASSERT(testShape(shapeB, Shape<Dna, GappedShape<HardwiredShape<1,1,3,1,2> > >(), true));

	pattern = "11110011";
	Shape<Dna, OneGappedShape> shapeC(pattern);
	SEQAN_TASSERT(testShape(shapeC, Shape<Dna, GappedShape<HardwiredShape<1,1,1,3,1> > >(), true));
}
