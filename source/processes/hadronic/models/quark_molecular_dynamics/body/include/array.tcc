#include <stdarg.h>
#include "g4std/iostream"


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
    array[i] = va_arg(ap,int);

  va_end(ap);
}

template<class t> 
Array<t>::Array(const Array<t>& x) : N(x.N),array(new t[x.N]) 
{
  for (int i=0; i<N; i++)
    array[i] = x[i];
}

