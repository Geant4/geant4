#include <stdarg.h>
#include <iostream.h>


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

template<class t> 
void Array<t>::insert(const t& x,int pos)
{
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
}
