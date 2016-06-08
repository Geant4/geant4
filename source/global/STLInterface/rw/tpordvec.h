#ifndef __tpordvec
#define __tpordvec

#ifndef G4USE_OLDSTL
# include <vector>
# include <set>
# include <functional>
#else
# include <vector.h>
# include <set.h>
# include <function.h>
#endif
#include "rw/defs.h"

#if defined(G4USE_NAMESPACE)
using std::vector;
using std::set;
#endif


template<class T> class RWTPtrOrderedVector : public vector<T*> {

typedef typename vector<T*>::iterator iterator;
typedef typename vector<T*>::const_iterator const_iterator;

public:
  RWTPtrOrderedVector(size_t capacity=RWDEFAULT_CAPACITY):
    vector<T*>(capacity,(T*)0),rwsize(0){}

  RWTPtrOrderedVector(const RWTPtrOrderedVector<T>& v):
    vector<T*>(v),rwsize(v.rwsize){}

  RWTPtrOrderedVector<T>& operator=(const RWTPtrOrderedVector<T>& v)
    {
      vector<T*>::operator=(v),
      rwsize=v.rwsize;
      return *this;
    }

  ~RWTPtrOrderedVector(){}

  T*& operator()(size_t i) 
    {
      //if(i>=rwsize)
      //	RWTHROW(RWBoundsErr("RWTPtrOrderedVector ()",rwsize,i));
      if(i<vector<T*>::size() && rwsize<=i) rwsize=i+1;
      return vector<T*>::operator[](i); 
    }

  T* const& operator()(size_t i) const 
    { 
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrOrderedVector const ()",rwsize,i));
      return vector<T*>::operator[](i); 
    }

  T*& operator[](size_t i) 
    {
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrOrderedVector []",rwsize,i));
      //if(i<vector<T*>::size() && rwsize<=i) rwsize=i+1;
      return vector<T*>::operator[](i); 
    }

  T* const& operator[](size_t i) const 
    {      
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrOrderedVector const []",rwsize,i));
       return vector<T*>::operator[](i); 
    }


  void clearAndDestroy();

  void clear();

  size_t index ( const T*) const;

  void resize(size_t);

  void insert ( T* a )
    { 
      if(rwsize<vector<T*>::size())
	{
	  vector<T*>::operator[](rwsize)=a;
	  rwsize++;
	}
      else
	{
	  vector<T*>::push_back(a);
	  rwsize=vector<T*>::size();
	}
    }

  //
  // insertAt should throw an exception of type RWBoundsErr if you try
  // to insert out of range
  //
  void insertAt ( size_t, T*); 

  size_t occurrencesOf(const T*) const;	

  T* remove ( const T*);	

  size_t removeAll ( const T* );	

  T* removeAt ( size_t);

  T* removeFirst ();

  T* removeLast ();

  T* find(const T*) const;															
  size_t entries() const { return rwsize; }

  void append ( T* a ) 
    { 
      insert( a );
    }

  size_t length() const { return rwsize; }

  T*& at ( size_t i ) 
    { 
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrOrderedVector at",rwsize,i));
      return operator()(i); 
    }

  T*const& at ( size_t i ) const 
    {       
      if(i>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrOrderedVector const at",rwsize,i));
      return operator()(i); 
    }

  RWBoolean contains ( const T* a) const;
  
  RWBoolean isEmpty()const {return rwsize==0;}
  
  T* first() const
    {
      if(!rwsize)
	RWTHROW(RWBoundsErr("RWTPtrOrderedVector front",rwsize,0));
      return vector<T*>::front();
    }

  T* last() const
    {
      if(rwsize)
	return vector<T*>::operator[](rwsize-1);
      else
	return 0;
    }

  void prepend(T* const ptr)
    {
      vector<T*>::insert(vector<T*>::begin(),ptr);
      rwsize++;
    }

private:

  size_t rwsize;


};


template<class T > void RWTPtrOrderedVector<T>::clear() 
{
  //vector<T*>::erase(vector<T*>::begin(), vector<T*>::end());
  rwsize=0;
}

template<class T > void RWTPtrOrderedVector<T>::clearAndDestroy() 
{
  int sz;
  set<T*,greater<T*> > tmp;
  for (sz=0;sz<rwsize;sz++)
    {
      T* current;
      current=vector<T*>::operator[](sz);
      if ( current )
         tmp.insert(current);
    }
  typename set<T*,greater<T*> >::iterator it;
  for(it=tmp.begin();it!=tmp.end();it++)
    {
      delete *it;
    }
  //vector<T*>::erase(vector<T*>::begin(), vector<T*>::end());
  rwsize=0;
}

template<class T> size_t RWTPtrOrderedVector<T>::index( const T* a ) const
{
  const_iterator i;
  size_t ptn = 0;
  for (i = vector<T*>::begin();ptn<rwsize; i++,ptn++)
    {
      if (**i==*a) return ptn;  
    }
  return (ptn=~(size_t)0);
}

template<class T>
size_t RWTPtrOrderedVector<T>::occurrencesOf( const T* a ) const  
{
  const_iterator i;
  size_t ptn = 0;
  size_t sz=0;
  for (i = vector<T*>::begin();sz<rwsize; i++,sz++)
    {
      if (**i==*a) ptn++;
    }
  return ptn;
}

template<class T> 
RWBoolean RWTPtrOrderedVector<T>::contains( const T* a ) const  
{
  const_iterator i;
  size_t ptn=0;
  for (i = vector<T*>::begin(); ptn<rwsize; i++,ptn++)
    {
      if (**i==*a) return true;
    }
  return false;
}

template<class T> T* RWTPtrOrderedVector<T>::remove( const T* a )  
{
  iterator i;
  size_t ptn=0;
  size_t rwsz=rwsize;
  for (i = vector<T*>::begin();ptn<rwsz; i++,ptn++)
    {
      if (**i==*a) 
	{
	  T* tmp=*i;
	  vector<T*>::erase(i);
	  rwsize--;
	  return tmp;
	} 
    }
  return 0;
}

template<class T> 
size_t RWTPtrOrderedVector<T>::removeAll( const T* a )  
{
  iterator i;
  size_t ptn = 0;
  size_t sz=0;
  size_t rwsz=rwsize;  
  for (i = vector<T*>::begin();sz<rwsz; i++,sz++)
    {
      if (**i==*a) 
	{
	  vector<T*>::erase(i);
	  rwsize--;
          i--;
	  ptn++;
	} 
    }
  return ptn;
}
// should throw exception if ...
template<class T> 
T* RWTPtrOrderedVector<T>::removeAt( size_t i )  
{
  if(i>=rwsize)
    RWTHROW(RWBoundsErr("RWTPtrOrderedVector removeAt",rwsize,i));
  iterator it=vector<T*>::begin();
  int j;
  for(j=0;j<rwsize && j<i;j++) it++;
  if(it!=vector<T*>::end())
    {
      T* tmp = operator[](i);
      vector<T*>::erase(it);
      rwsize--;
      return tmp;
    }
  else
    return 0;
}

template<class T>	
T* RWTPtrOrderedVector<T>::find(const T* a) const 
{
  const_iterator i;
  size_t ptn=0;
  for (i = vector<T*>::begin(); ptn<rwsize; i++,ptn++)
    {
      if (**i==*a) 
	{
	  return *i;
	} 
    }
  return 0;
}															
template<class T>	
void RWTPtrOrderedVector<T>::resize(size_t N)
{
  if(N>vector<T*>::size())
    {
      int e=N-vector<T*>::size();
      for(int i=0;i<e;i++)
	vector<T*>::push_back(0);
    }
  else if(N>rwsize)
    {
      iterator it=vector<T*>::begin();
      for(int i=0;i<N;i++,it++);
      vector<T*>::erase(it,vector<T*>::end());
    }
}

template<class T>	
void RWTPtrOrderedVector<T>::insertAt ( size_t i, T* a ) 
{
  if(i>rwsize)
    RWTHROW(RWBoundsErr("RWTPtrOrderedVector insertAt",rwsize,i));
  if(rwsize<vector<T*>::size())
    {
      for(int j=rwsize;j>i;j--)
	vector<T*>::operator[](j)=vector<T*>::operator[](j-1);
      vector<T*>::operator[](i)=a;
      rwsize++;
    }
  else
    {
      iterator it=vector<T*>::begin();
      for(int j=0;j<i;j++,it++);
      vector<T*>::insert(it,a);
      rwsize++;
    }
}

template<class T>	
T* RWTPtrOrderedVector<T>::removeFirst ()
{
  if(rwsize!=0)
    {
      T* tmp=vector<T*>::front();
      vector<T*>::erase(vector<T*>::begin());
      rwsize--;
      return tmp;
    }
  else
    return 0;
}

template<class T>
T* RWTPtrOrderedVector<T>::removeLast ()
{
  if(rwsize==0) return (T*) 0;
  
  rwsize--;
  T* tmp=vector<T*>::operator[](rwsize);
  vector<T*>::operator[](rwsize)=0;

  return tmp;
}


#endif

