#ifndef __tphdict
#define __tphdict

#include "rw/defs.h"
#ifndef G4USE_OLDSTL
  #include <vector>
#else
  #include <vector.h>
#endif

#ifndef G4USE_OLDSTL
#include <hash_map>
#include <set>
#include <functional>
#else
# include <set.h>
# include <function.h>
# if defined(WIN32) && defined(OS_WIN_NT_4_0)
//   OBJECT Space on WIN32 choose to change this name..
#   include <hashmap.h>
# else
#   include <hash_map.h>
# endif
#endif

template <class PK>
class Equal_Ptr_To
{
public:
  Equal_Ptr_To(){}
  RWBoolean operator()(const PK*a,const PK*b) const
    {return *a==*b;}
};
  
#if defined(G4USE_NAMESPACE)
using std::set;
#endif

#include "rw/rwstlhash.h"


template < class PK >
class HashPtr {
public:
  HashPtr(unsigned (*f)( const PK &)=HashDefault):fhashfun(f){}
  HashPtr(const HashPtr<PK>&s):fhashfun(s.fhashfun){}
  unsigned operator()( PK* key ) const;
  void SetHashFun( unsigned (*hashfun)( const PK& ) );
private:
  unsigned (*fhashfun)( const PK & );  
};

#define hash_dict_t hash_map< K*, V*, HashPtr<K>,Equal_Ptr_To<K> >

template < class K, class V >
class RWTPtrHashDictionary 
{

#ifdef G4USE_EXPLICIT_TYPES_IN_TEMPLATES
  typedef os_hash_table_iterator
    <
    OS_PAIR( K*, V* ),
    OS_DIFF_TYPE( OS_PAIR( K*, V* ) )
    >  iterator;
  typedef os_hash_table_const_iterator
    <
    OS_PAIR( K*, V* ),
    OS_DIFF_TYPE( OS_PAIR( K*, V* ) )
    >  const_iterator;
#else
  typedef hash_dict_t::iterator iterator;
  typedef hash_dict_t::const_iterator const_iterator;
#endif

public:

  RWTPtrHashDictionary( unsigned (*hashfun)( const K& key ),size_t n=100)
    :fhm(n,HashPtr<K>(hashfun),Equal_Ptr_To<K>()){}
    
  //  Default copy constructor ?
  RWTPtrHashDictionary( const RWTPtrHashDictionary<K,V>& aPtrHashDict ):
    fhm( aPtrHashDict.fhm ){} 

  V*& operator[]( K*& key ) {return fhm[key]; }
  void clear() { fhm.clear(); }
  void clearAndDestroy( void ) 
    {
      //We are deleting in reverse order
      set<K*,greater<K*> > ktmp;
      set<V*,greater<V*> > vtmp;
      for(iterator it=fhm.begin();
	  it!=fhm.end();it++)
	{
	  K* kcurrent;
	  V* vcurrent;
	  kcurrent=((*it).first);
	  vcurrent=((*it).second);
	  if ( kcurrent )
	    ktmp.insert(kcurrent);
	  if ( vcurrent )
	    vtmp.insert(vcurrent);
	}
      typename set<K*,greater<K*> >::iterator kit;
      for(kit=ktmp.begin();kit!=ktmp.end();kit++)
	{
	  delete *kit;
	}
      typename set<V*,greater<V*> >::iterator vit;
      for(vit=vtmp.begin();vit!=vtmp.end();vit++)
	{
	  delete *vit;
	}
      fhm.clear();
    }

  void clearAndDestroyKey(void)
    {
      //We are deleting in reverse order
      set<K*,greater<K*> > ktmp;
      for(iterator it=fhm.begin();
	  it!=fhm.end();it++)
	{
	  K* kcurrent;
	  kcurrent=((*it).first);
	  if ( kcurrent )
	    ktmp.insert(kcurrent);
	}
      typename set<K*,greater<K*> >::iterator kit;
      for(kit=ktmp.begin();kit!=ktmp.end();kit++)
	{
	  delete *kit;
	}
      fhm.clear();
    }

  K* find ( const K* key )
    {
      iterator it;
      // The STL RW const problem
      K tval(*key);
      K* tptr=&tval;
      it=fhm.find(tptr);
      if(it!=fhm.end())
        return (*it).first;
      else
        return 0;
    }

  V* findValue( const K* key ) 
    { 
      iterator it;
      // The STL RW const problem
      K tval(*key);
      K* tptr=&tval;
      it=fhm.find(tptr);
      if(it!=fhm.end())
	return (*it).second;
      else
	return 0;
    }

  void insertKeyAndValue (K* key, V* value) {fhm[key]=value;}
  K* remove(const K* key) { 
    iterator i;
    K tmp(*key);
    K* ptr=&tmp;
    i = fhm.find(ptr);
    if ( i == fhm.end() ) { return NULL; }
    else {
      ptr=(*i).first;
      fhm.erase( i );
      return ptr;
    }    
  };
  void resize( unsigned N ) { fhm.resize(N); } 
  K* findKeyAndValue( const K* key, V*& retval ) 
    { 
      iterator i;
      K tmp(*key);
      K* ptr=&tmp;
      i = fhm.find(ptr );
      if ( i == fhm.end() ) { return false; }
      else {
	retval = (*i).second; 
	return (*i).first;
      }
    }

  size_t entries() const {return fhm.size();}
  size_t max_size() const {return fhm.max_size();}
  size_t bucket_count() const {return fhm.bucket_count();}
  
  RWBoolean contains(const K* Key) const
    {
      const_iterator it;
      // The STL RW const problem
      K tval(*Key);
      K* tptr=&tval;
      it=fhm.find(tptr);
      if(it!=fhm.end())
	return true;
      else
	return false;
    }


  //This is STL, so will not be used by RW 
  iterator begin(){return fhm.begin();}
  iterator end(){return fhm.end();}

private:

  hash_dict_t fhm;
};

template < class K, class V >
class RWTPtrHashDictionaryIterator {

#ifdef G4USE_EXPLICIT_TYPES_IN_TEMPLATES
  typedef os_hash_table_iterator
    <
    OS_PAIR( K*, V* ),
    OS_DIFF_TYPE( OS_PAIR( K*, V* ) )
    >  iterator;
  typedef os_hash_table_const_iterator
    <
    OS_PAIR( K*, V* ),
    OS_DIFF_TYPE( OS_PAIR( K*, V* ) )
    >  const_iterator;
#else
  typedef hash_dict_t::iterator iterator;
  typedef hash_dict_t::const_iterator const_iterator;
#endif

public:
  RWTPtrHashDictionaryIterator(RWTPtrHashDictionary<K,V>&adict):
    mydict(&adict),it(adict.begin()),defined(false){}

  RWBoolean operator++ ()
    {
      if(!defined) return false;
      it++;
      return it!=mydict->end() ? true : false;
    }
      
  RWBoolean operator()()
    {
      if(defined)
	return operator++();
      else
	{
	  defined=true;
	  it=mydict->begin();
	  return it!=mydict->end() ? true : false;
	}
    }

  K* key () const {return (*it).first;}

  void reset (){defined=false;}
  V* value () const{return (*it).second;}

private:


  iterator it;  
  RWTPtrHashDictionary<K,V>* mydict;
  RWBoolean defined;
};

template< class PK >
unsigned HashPtr<PK>::operator()( PK * key ) const {
    return fhashfun( *key );
}

template< class PK > 
void HashPtr<PK>::SetHashFun( unsigned (*hashfun)( const PK& ) ) {
  fhashfun = hashfun;
}

#undef hash_dict_t

#endif













