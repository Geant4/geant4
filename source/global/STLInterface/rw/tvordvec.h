#ifndef __tvordvec
#define __tvordvec

#ifndef G4USE_OLDSTL
  #include <vector>
  #include <algorithm>
#else
  #include <vector.h>
  #include <algo.h>
#endif

#include "rw/defs.h"
#include "rw/tvvector.h"

template<class T> class RWTValOrderedVector : public vector<T> 
{
 
typedef typename vector<T>::iterator iterator;
typedef typename vector<T>::const_iterator const_iterator;
 
public:
  RWTValOrderedVector(size_t capacity=RWDEFAULT_CAPACITY):vector<T>(capacity),rwsize(0){}

  RWTValOrderedVector(const RWTValOrderedVector<T>& v):
    vector<T>(v),rwsize(v.rwsize){}


  RWTValOrderedVector<T>& operator=(const RWTValOrderedVector<T>& v)
    {
      vector<T>::operator=(v);
      rwsize=v.rwsize;
      return *this;
    }

 ~RWTValOrderedVector(){}

  T& operator()(size_t i) 
    {
      //if(i>=rwsize)
      //	RWTHROW(RWBoundsErr("RWTValOrderedVector ()",rwsize,i));
      if(i<vector<T>::size() && rwsize<=i) rwsize=i+1;
      return vector<T>::operator[](i); 
    }

  const T& operator()(size_t i) const
   { 
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTValOrderedVector const()",rwsize,i));
      return vector<T>::operator[](i);
   }

  T& operator[](size_t i) 
    {
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTValOrderedVector []",rwsize,i));
//      if(i<vector<T>::size() && rwsize<=i) rwsize=i+1;
      return vector<T>::operator[](i); 
    }

  const T& operator[](size_t i) const 
    { 
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTValOrderedVector const []",rwsize,i));
      return vector<T>::operator[](i); 
    }

  T& at(size_t n )
    {
      if(n>=rwsize)
	RWTHROW(RWBoundsErr("RWTValOrderedVector at",rwsize,n));
      return operator()(n);
    }

  T  at(size_t n ) const
    {
      if(n>=rwsize)
	RWTHROW(RWBoundsErr("RWTValOrderedVector const at",rwsize,n));
      return operator()(n);
    }

  void resize(size_t);

  void clear( void ) 
    { 
      //vector<T>::erase(vector<T>::begin(), vector<T>::end());
      for(size_t i=0;i<rwsize;i++)
      vector<T>::operator[](i)=T();
      rwsize=0;
    }

  RWBoolean contains( const T) const;  

  size_t length() const{ return rwsize;}

  size_t index ( const T&);

  RWBoolean insert ( const T&); 

  //
  // insertAt should throw an exception of type RWBoundsErr if you try
  // to insert out of range
  //
  void insertAt ( size_t, const T&); 

  size_t occurrencesOf ( const T& a ) const;	

  RWBoolean remove ( const T& a );	

  size_t removeAll ( const T& a );	

  T removeAt ( size_t i );
	
  size_t entries () const { return rwsize; }

  void append ( const T& a )
    { 
      insert ( a );
    }

  T first() const
    {
      if(!rwsize)
	RWTHROW(RWBoundsErr("RWTValOrderedVector first",rwsize,0));
      return vector<T>::front();
    }

  T last() const
    {
      if(rwsize)
	return vector<T>::operator[](rwsize-1);
      else
	RWTHROW(RWBoundsErr("RWTValOrderedVector last",rwsize,0));
    }

private:

  size_t rwsize;
  
};

template<class T> 
size_t RWTValOrderedVector<T>::index( const T& a ) 
{
  iterator i;
  size_t ptn = 0;
  for (i = vector<T>::begin();ptn<rwsize; i++,ptn++)
    {
      if (*i==a) return ptn;  
    }
  return (ptn=~(size_t)0);
}

template<class T> 
size_t RWTValOrderedVector<T>::occurrencesOf( const T& a ) const  
{
  const_iterator i;
  size_t ptn = 0;
  size_t sz=0;
  for (i = vector<T>::begin(); sz<rwsize; i++,sz++)
    {
      if (*i==a) ptn++;
    }
  return ptn;
}

template<class T> 
RWBoolean RWTValOrderedVector<T>::remove( const T& a )  
{
  iterator i;
  size_t sz=0;
  size_t rwsz=rwsize;
  for (i = vector<T>::begin(); sz<rwsz; i++,sz++)
    {
      if (*i==a) 
	{
	  vector<T>::erase(i);
	  rwsize--;
	  return true;
	}
    } 
  return false;
}

template<class T> size_t RWTValOrderedVector<T>::removeAll
( const T& a )  
{
  iterator i;
  size_t ptn = 0;
  size_t sz=0;
  size_t rwsz=rwsize;
  for (i = vector<T>::begin();sz<rwsz; i++,sz++)
    {
      if (*i==a) 
	{
	  vector<T>::erase(i);
	  ptn++;
          i--;
	  rwsize--;
	} 
    }
  return ptn;
}

template<class T>
T RWTValOrderedVector<T>::removeAt( size_t i )  
{
  if(i>=rwsize)
    RWTHROW(RWBoundsErr("RWTValOrderedVector removeAt",rwsize,i));
  T tmp = operator[](i);
  iterator it=vector<T>::begin();
  int j;
  for(j=0;j<rwsize && j<i;j++) it++;
  if(it!=vector<T>::end())
    {
      vector<T>::erase(it);
      rwsize--;
    }
  return tmp;
}



template<class T>
void RWTValOrderedVector<T>::resize(size_t N)
{
  if(N>vector<T>::size())
    {
      int e=N-vector<T>::size();
      for(int i=0;i<e;i++)
	vector<T>::push_back(T());
    }
  else if(N>rwsize )
    {
      iterator it=vector<T>::begin();
      for(int i=0;i<N;i++,it++);
      vector<T>::erase(it,vector<T>::end());
    }
}


template<class T>
RWBoolean RWTValOrderedVector<T>::contains( const T a ) const  
{
  const_iterator i;
  size_t sz=0;
  for (i = vector<T>::begin();sz<rwsize; i++,sz++)
    {
      if (*i==a) return true;
    }
  return false;
} 

template<class T>
RWBoolean RWTValOrderedVector<T>::insert ( const T& a ) 
{ 
  if(rwsize<vector<T>::size())
    {
      this->vector<T>::operator[](rwsize)=a;
      rwsize++;
    }
  else
    {
      vector<T>::push_back(a);
      rwsize=vector<T>::size();
    }
  return true;
}

template<class T>
void RWTValOrderedVector<T>::insertAt ( size_t i, const T& a ) 
{
  if(i>rwsize)
    RWTHROW(RWBoundsErr("RWTValOrderedVector insertAt",rwsize,i));
  if(rwsize<vector<T>::size())
    {
      for(int j=rwsize;j>i;j--)
	vector<T>::operator[](j)=this->vector<T>::operator[](j-1);
      vector<T>::operator[](i)=a;
      rwsize++;
    }
  else
    {
      iterator it=vector<T>::begin();
      for(int j=0;j<i;j++,it++);
      vector<T>::insert(it,a);
      rwsize++;
    }
}
#endif

