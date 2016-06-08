#ifndef __tpvector
#define __tpvector

#ifndef G4USE_OLDSTL
#include <vector>
#else
#include <vector.h>
#endif

#if defined(G4USE_NAMESPACE)
using std::vector;
#endif

template<class T> class RWTPtrVector : public vector<T*> {
 
typedef typename vector<T*>::iterator iterator;
typedef typename vector<T*>::const_iterator const_iterator;
 
public:

  RWTPtrVector ();
  RWTPtrVector (unsigned int n);
  RWTPtrVector (unsigned int n, T* const & iPtr);
  RWTPtrVector (const RWTPtrVector<T>& iPtr);
  ~RWTPtrVector(){}
  // Copy constructor and assignment operator inherited from vector<T*>.
  RWTPtrVector<T>& operator = (T* p) {
    for(iterator i = vector<T*>::begin();
	i != vector<T*>::end();
	++i) {
	  *i = p;
	}
    return *this;
  }
	T* operator () ( size_t n ) const { return vector<T*>::operator[](n);}
	T*& operator () ( size_t n ) { return vector<T*>::operator[](n);}
//	
// The [] operator should perform bound checking
//	
 	T* operator [] ( size_t n ) const { return vector<T*>::operator[](n);}
 	T*& operator [] ( size_t n ) { return vector<T*>::operator[](n);}
	T* const * data() const { return &(vector<T*>::front()); }
	size_t length() { return vector<T*>::size();}
	void reshape( size_t n ); 
};

#include "rw/tpvector.icc"

#endif
