#ifndef __array__
#define __array__

#ifdef IS_GCC
#pragma interface
#endif

#include <iostream.h>

template<class t>
class Array
{
  friend ostream& operator<<(ostream& o,const Array<t>& i){
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
  Array();
  Array(int n);
  Array(int n,t* x);
  Array(int n,t a1 ...);
  Array(const Array<t>&);
  ~Array() { delete [] array; }
  operator t*() { return array; }
  t operator[](int i) const { return array[i]; }
  t& operator[](int i) { return array[i]; }
  void insert(const t& x,int pos = -1);
  int size() const { return N; }
};

#ifndef IS_GCC
#include "array.tcc"
#endif

#endif
