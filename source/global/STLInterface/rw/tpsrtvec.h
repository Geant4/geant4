#ifndef __tpsrtvec
#define __tpsrtvec

#include "rw/rw2stl_algo.h"

#ifndef G4USE_OLDSTL
  #include <vector>
  #include <algorithm>
#else
  #include <vector.h>
  #include <algo.h>
#endif

#if defined(G4USE_NAMESPACE)
using std::vector;
using std::sort;
#endif

#include "rw/defs.h"


template <class T>
class RWTPtrSortedVector : public vector<T*> {
 
#ifdef G4USE_EXPLICIT_TYPES_IN_TEMPLATES
  typedef T** iterator;
  typedef T* const* const_iterator;
#endif
 
public:
  RWTPtrSortedVector(size_t n=RWDEFAULT_CAPACITY):vector<T*>(n),rwsize(0)
    {
      for(int i=0;i<n;i++)
	vector<T*>::operator[](i)=0;
    }

  RWTPtrSortedVector(const RWTPtrSortedVector<T>& v):vector<T*>(v),
    rwsize(v.rwsize){}

  RWTPtrSortedVector<T>& operator=(const RWTPtrSortedVector<T>& v)
  {
      vector<T*>::operator=(v);
      rwsize=v.rwsize;
      return *this;     
  }
  

  ~RWTPtrSortedVector(){}

  T* const& operator [] (size_t n) const
    {
      if(n<0 || n>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrSortedVector",rwsize,n));
      return vector<T*>::operator[](n);
    }
  T* const& operator () (size_t n) const
    {
      return vector<T*>::operator[](n);
    }

  T*& at(size_t n )
    {
      if(n<0 || n>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrSortedVector",rwsize,n));
      if(n<vector<T*>::size() && rwsize<=n) rwsize=n+1;
      return vector<T*>::operator[](n);
    }

  T*  at(size_t n ) const
    {
      if(n<0 || n>=rwsize)
	RWTHROW(RWBoundsErr("RWTPtrSortedVector",rwsize,n));
      return vector<T*>::operator[](n);
    }

  void clear()
    {
      vector<T*>::erase(vector<T*>::begin(),vector<T*>::end());
      rwsize=0;
    }

  void clearAndDestroy () 
    { 
      size_t sz=0;
      for(iterator it=vector<T*>::begin();sz<rwsize;it++,sz++)
	delete *it;
      vector<T*>::erase(vector<T*>::begin(),vector<T*>::end());
      rwsize=0;
    }

  T* find (const T* a) const {
    if(!a || vector<T*>::empty()) return 0;
    //We are sorted, so try a n log n search
    int t,d;
    int u=rwsize,l=0;
    do
      {
	d=(u-l)/2;
	t=d+l;
	//this asks *a <= (*(*this)[t] whithout requiring operator<=
	if(sorter(a,(*this)[t]) || *a==*(*this)[t])
	  u=t;
	else
	  l=t;
      }
    while(d>0);
    
    if(*(*this)[u]==*a)
      return (*this)[u];
    else
      return 0;
  }

  size_t entries () const {return rwsize;}

  void insert ( T* a ) 
    { 
      if(rwsize<vector<T*>::size())
	{
	  //We have empty entries
	  iterator it=vector<T*>::begin();
	  for(int i=0;i<rwsize;i++,it++);
	  *it=a;
	  it++;
	  sort(vector<T*>::begin(),it,sorter);
	}
      else
	{
	  vector<T*>::push_back(a);
	  //This makes it a sorted thing;
	  sort(vector<T*>::begin(),vector<T*>::end(),sorter);
	}
      rwsize++;
    }

  size_t index ( const T* a ) const 
    {
      T*const* ii = vector<T*>::begin();
      size_t cnt = 0;
      for (T*const* it = ii;cnt<rwsize; ++it,cnt++) 
	{
	  if ( **it == *a ) return cnt;
	}
      return 0;
    }

  RWBoolean isEmpty () const{return rwsize==0;}

  T* const & first () const 
    {
      if(!rwsize)
	RWTHROW(RWBoundsErr("RWTPtrSortedVector first()",rwsize,0));
      return vector<T*>::front();
    }
  
  T* const & last () const 
    {
      if(rwsize)
	return vector<T*>::operator[](rwsize-1);
      else
	RWTHROW(RWBoundsErr("RWTPtrSortedVector last()",rwsize,0));
    }

  size_t occurrencesOf(const T* a) const 
    {
      size_t cnt = 0;
      size_t sz = 0;
      for (const_iterator it = vector<T*>::begin();sz<rwsize; ++it,sz++) 
	{
	  if ( **it == *a ) cnt++;
	}
      return cnt;
    }

  T* remove ( const T* a )
    {
      iterator retval;
      size_t sz=0;
      for(retval=vector<T*>::begin();sz<rwsize;retval++,sz++)
	if(**retval==*a) break;
      if(sz<rwsize && retval!=vector<T*>::end()) 
	{
	  T* tmp=*retval;
	  vector<T*>::erase(retval);
	  rwsize--;
	  return tmp;
	}
      else
	return 0;
    }

  size_t removeAll ( const T* a)
    {
      iterator s=vector<T*>::end();
      iterator e=vector<T*>::end();
      iterator r;
      size_t sz=0;
      int i=0;
      for(r=vector<T*>::begin();sz<rwsize;r++,sz++)
	{
	  if(i==0)
	    {
	      if(**r==*a)
		{
		  s=r;
		  i++;
		}
	    }
	  else 
	    {
	      if(**r!=*a)
		{
		  e=r;
		  break;
		}
	      else
		i++;
	    }
	}
      vector<T*>::erase(s,e);
      rwsize-=i;
      return i;
    }


private:

  size_t rwsize; 
  RW2STL_LessPtr<T> sorter;  
	
};
#endif









