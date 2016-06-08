#ifndef __tvvector
#define __tvvector

#ifndef G4USE_OLDSTL
  #include <vector>
#else
  #include <vector.h>
#endif

#if defined(G4USE_NAMESPACE)
using std::vector;
#endif

#include "rw/defs.h"

template<class T> class RWTValVector : public vector<T> {
 
typedef typename vector<T>::iterator iterator;
typedef typename vector<T>::const_iterator const_iterator;
 
public:
  RWTValVector ():rwsize(0){}

  RWTValVector (size_t n):vector<T>(n),rwsize(n){}

  RWTValVector (size_t n, const T& ival):vector<T>(n,ival),rwsize(n){}

  RWTValVector (const RWTValVector<T>& ival):vector<T>(ival),rwsize(ival.rwsize){}

  ~RWTValVector(){};

  RWTValVector<T>& operator= ( const RWTValVector<T>& v)
    {
      vector<T>::operator=(v);
      rwsize=v.rwsize;
      return *this;     
    }

  size_t length() const { return vector<T>::size();}

  T& operator()(size_t i) 
    {
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTValVector operator()",rwsize,i));
      if(i<vector<T>::size() && rwsize<i) rwsize=i+1;
      return vector<T>::operator[](i); 
    }

  const T& operator()(size_t i) const { return vector<T>::operator[](i); }

  T& operator[](size_t i) 
    {
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTValVector operator[]",rwsize,i));
      if(i<vector<T>::size() && rwsize<i) rwsize=i+1;
      return vector<T>::operator[](i); 
    }

  const T& operator[](size_t i) const 
    { 
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTValVector const operator[]",rwsize,i));
      return vector<T>::operator[](i); 
    }

  const T* data() const { return &(vector<T>::front()); };

  void resize(size_t N)
    {
      if(N>vector<T>::size())
	{
	  int e=N-vector<T>::size();
	  for(int i=0;i<e;i++)
	    vector<T>::push_back(T());
	}
      else
	{
	  iterator it=vector<T>::begin();
	  for(int i=1;i<=N;i++,it++);
	  vector<T>::erase(it,vector<T>::end());
	}
       rwsize=N;
    }

  void reshape( size_t n )
    {
      resize(n);
    }

  const T& ref(size_t i) const{return vector<T>::operator[](i);}	

private:

  size_t rwsize;

};
#endif

