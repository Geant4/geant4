#ifndef __array__
#define __array__

#ifdef IS_GCC
#pragma interface
#endif

#include <stdarg.h>
#include "g4std/iostream"


template<class t> class Array;

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
  Array();
  Array(int n);
  Array(int n,t* x);
  Array(int n,t a1 ...);
  Array(const Array<t>&);
  ~Array() { delete [] array; }
  operator t*() { return array; }
  t operator[](int i) const { return array[i]; }
  t& operator[](int i) { return array[i]; }
  void insert(const t& x,int pos = -1) {
// ------------------------------------------
// implementation from array.tcc:
    if ( pos<0 ) pos = N;
    t* new_array = new t[N+1];
    for (int i=0; i<pos; i++)
      new_array[i] = array[i];
    new_array[pos] = x;
    for (int j=pos; j<N; j++)
      new_array[j+1] = array[j];
    if ( array ) 
      delete [] array;
    array = new_array;
    ++N;
// ------------------------------------------
  }
  int size() const { return N; }
};

// ------------------------------------------
// implementation from array.tcc:

template<class t> 
Array<t>::Array() : N(0),array(0) {}

template<class t> 
Array<t>::Array(int n) : N(n),array(new t[n]) {}

template<class t> 
Array<t>::Array(int n,t* x) : N(n),array(new t[n]) 
{
  for (int i=0; i<N; i++)
    array[i] = x[i];
}

template<class t> 
Array<t>::Array(int n,t a1 ...) : N(n),array(new t[n]) 
{
  va_list ap;
  va_start(ap,a1);

  array[0] = a1;
  for (int i=1; i<N; i++)
    array[i] = va_arg(ap,t);

  va_end(ap);
}

template<class t> 
Array<t>::Array(const Array<t>& x) : N(x.N),array(new t[x.N]) 
{
  for (int i=0; i<N; i++)
    array[i] = x[i];
}

// ------------------------------------------

#endif
