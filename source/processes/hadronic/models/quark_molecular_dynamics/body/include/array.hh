#ifndef __array__
#define __array__

#ifdef IS_GCC
#pragma interface
#endif

#include "g4std/iostream"

template<class t>
class Array
{
  friend G4std::ostream& operator<<(G4std::ostream& o,const Array<t>& i){
    o << "(";
    for (int k=0; k<i.N; k++) {
      o << i.array[k];
      if ( k<i.N-1 ) 
	o << ",";
    }
    o << ")";
    return o;
  }
  int N;
  t* array;
public:
  Array(int n);
  Array(int n,t* x);
  Array(int n,t a1 ...);
  Array(const Array<t>&);
  ~Array() { delete [] array; }
  operator t*() { return array; }
  t operator[](int i) const { return array[i]; }
  t& operator[](int i) { return array[i]; }
};

#ifndef IS_GCC
#include "array.tcc"
#endif

#endif
