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

#ifndef SEQAN_HEADER_RNA_ALPHABET_H
#define SEQAN_HEADER_RNA_ALPHABET_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// RNA5 Alphabet
//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Rna5_2_Ascii
{
	static char const VALUE[5];
};
template <typename T>
char const _Translate_Table_Rna5_2_Ascii<T>::VALUE[5] = {'A', 'C', 'G', 'U', 'N'}; 

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Rna_2_Ascii
{
	static char const VALUE[4];
};
template <typename T>
char const _Translate_Table_Rna_2_Ascii<T>::VALUE[4] = {'A', 'C', 'G', 'U'}; 


//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Byte_2_Rna5
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_Rna5<T>::VALUE[256] = 
{
	0,   1,   2,   3,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_Byte_2_Rna
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Byte_2_Rna<T>::VALUE[256] = 
{
	0,   1,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};


//////////////////////////////////////////////////////////////////////////////


template <typename T = void>
struct _Translate_Table_Ascii_2_Rna5
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_Rna5<T>::VALUE[256] = 
{
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //0
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //1
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //2
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //3

	4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	4,   4,   4,   4,   4,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	4,   0,   4,   1,   4,   4,   4,   2,   4,   4,   4,   4,   4,   4,   4,   4, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	4,   4,   4,   4,   4,   3,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //8
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //9
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //10
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //11
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //12
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //13
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4, //14
	4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4,   4  //15
};

//////////////////////////////////////////////////////////////////////////////


template <typename T = void>
struct _Translate_Table_Ascii_2_Rna
{
	static char const VALUE[256];
};
template <typename T>
char const _Translate_Table_Ascii_2_Rna<T>::VALUE[256] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   1,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0, //6
//   ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

struct _Rna {};
typedef SimpleType<unsigned char,_Rna> Rna;

template <> struct ValueSize< Rna > { enum { VALUE = 4 }; };
template <> struct BitsPerValue< Rna > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////

struct _Rna5 {};
typedef SimpleType<unsigned char, _Rna5> Rna5;

template <> struct ValueSize< Rna5 > { enum { VALUE = 5 }; };
template <> struct BitsPerValue< Rna5 > { enum { VALUE = 3 }; };

//////////////////////////////////////////////////////////////////////////////
//Rna assignment
//////////////////////////////////////////////////////////////////////////////

inline void 
assign(Ascii& target,
	   Rna const & source)
{
	SEQAN_CHECKPOINT
	target = _Translate_Table_Rna_2_Ascii<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<Rna, Byte> { typedef Rna Type; };
inline void assign(Rna & target, Byte c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_Rna<>::VALUE[c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<Rna, Ascii> { typedef Rna Type; };
inline void assign(Rna & target, Ascii c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Rna<>::VALUE[(unsigned char)c_source];
}

//////////////////////////////////////////////////////////////////////////////


template <>
struct CompareType<Rna, Unicode> { typedef Rna Type; };
inline void assign(Rna & target, Unicode c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Rna<>::VALUE[(unsigned char) c_source];
}

//____________________________________________________________________________

template <>
struct CompareType<Rna, Rna5> { typedef Rna Type; };
inline void assign(Rna & target, Rna5 const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value & 0x03;
}

//////////////////////////////////////////////////////////////////////////////
//Rna5 assignment
//////////////////////////////////////////////////////////////////////////////

inline void 
assign(Ascii& target,
	   Rna5 const & source)
{
	SEQAN_CHECKPOINT
	target = _Translate_Table_Rna5_2_Ascii<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<Rna5, Byte> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Byte c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_Rna5<>::VALUE[c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<Rna5, Ascii> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Ascii c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Rna5<>::VALUE[(unsigned char)c_source];
}

//////////////////////////////////////////////////////////////////////////////


template <>
struct CompareType<Rna5, Unicode> { typedef Rna5 Type; };
inline void assign(Rna5 & target, Unicode c_source)
{
	SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_Rna5<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////

template <>
struct CompareType<Rna5, Rna> { typedef Dna Type; };
inline void assign(Rna5 & target, Rna const & c_source)
{
SEQAN_CHECKPOINT
	target.value = c_source.value;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
