#ifndef __tvhdict
#define __tvhdict

#include "rw/defs.h"

#ifndef G4USE_OLDSTL
#include <hash_map>
#else
#  if defined(WIN32) && defined(OS_WIN_NT_4_0)
//   OBJECT Space on WIN32 choose to change this name..
#    include <hashmap.h>
#  else
#    include <hash_map.h>
#  endif
#endif


#include "rw/rwstlhash.h"

template < class K >
class Hash {
public:
  Hash(unsigned (*f)( const K& )=HashDefault):fhashfun(f){}
  unsigned operator()( const K& key ) const;
  void SetHashFun( unsigned (*hashfun )( const K& ) );
private:
  unsigned (*fhashfun)( const K& );  
};

#define hash_dict_t  hash_map< K, V, Hash<K>,equal_to<K> >

template < class K, class V >
class RWTValHashDictionary {

#ifdef G4USE_EXPLICIT_TYPES_IN_TEMPLATES
  typedef os_hash_table_iterator
    <
    OS_PAIR( K, V ),
    OS_DIFF_TYPE( OS_PAIR( K, V ) )
    >  iterator;
  typedef os_hash_table_const_iterator
    <
    OS_PAIR( K, V ),
    OS_DIFF_TYPE( OS_PAIR( K, V ) )
    >  const_iterator;
#else
  typedef hash_dict_t::iterator iterator;
  typedef hash_dict_t::const_iterator const_iterator;
#endif

public:

  RWTValHashDictionary(){}

  RWTValHashDictionary( unsigned (*hashfun)( const K& key) ,size_t n=100):
    fhm(n,Hash<K>(hashfun),equal_to<K>()){}

  //pair<iterator, RWBoolean> insert(const pair<K,V>& x);
  V& operator[]( const K& key ) { return fhm[key]; }
  void clear() { fhm.clear(); }
  RWBoolean findValue( const K& key, V& retval ) 
    { 
      iterator i;
      i = fhm.find( key );
      if ( i == fhm.end() ) { return false; }
      else {
	retval = (*i).second; 
	return true;
      }
    }

  RWBoolean contains(const K& key) const 
    {
      const_iterator i=fhm.find(key);
      if(i==fhm.end()) return false;
      return true;
    }

  RWBoolean remove(const K& key)
    {
      iterator i=fhm.find(key);
      if(i==fhm.end()) return false;
      fhm.erase(i);
      return true;
    }

  RWBoolean isEmpty() const
    {
      if(fhm.empty()) return true;
      return false;
    }

  size_t entries()const{return fhm.size();}

  void insertKeyAndValue(const K&key,const V&value)
    {
      fhm[key]=value;
    }

  //This is STL, so will not be used by RW 
  iterator begin(){return fhm.begin();}
  iterator end(){return fhm.end();}


private:
  hash_dict_t fhm;
};

template < class K, class V >
class RWTValHashDictionaryIterator {

#ifdef G4USE_EXPLICIT_TYPES_IN_TEMPLATES
  typedef os_hash_table_iterator
    <
    OS_PAIR( K, V ),
    OS_DIFF_TYPE( OS_PAIR( K, V ) )
    >  iterator;
  typedef os_hash_table_const_iterator
    <
    OS_PAIR( K, V ),
    OS_DIFF_TYPE( OS_PAIR( K, V ) )
    >  const_iterator;
#else
  typedef hash_dict_t::iterator iterator;
  typedef hash_dict_t::const_iterator const_iterator;
#endif
 
public:
  RWTValHashDictionaryIterator(RWTValHashDictionary<K,V>&adict):
    mydict(&adict),it(adict.begin()),defined(false){}

  RWBoolean operator++()
    {
      if(!defined) return false;
      it++;
      return it!=mydict->end() ? true : false;
    }

  RWBoolean operator++(int)
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
      
  K key () const {return (*it).first;}

  void reset (){defined=false;}
  void reset (RWTValHashDictionary<K,V>&adict)
    {
      mydict=adict;
      defined=false;
    }
  V value () const{return (*it).second;}

private:
  iterator it;  
  RWTValHashDictionary<K,V>* mydict;
  RWBoolean defined;
};

#include "rw/tvhdict.icc"

#undef hash_dict_t

#endif
