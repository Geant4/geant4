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
class G4Hash {
public:
  G4Hash(unsigned (*f)( const K& )=HashDefault):fhashfun(f){}
  unsigned operator()( const K& key ) const;
  void SetHashFun( unsigned (*hashfun )( const K& ) );
private:
  unsigned (*fhashfun)( const K& );  
};

#define hash_dict_t  hash_map< K, V, G4Hash<K>,equal_to<K> >

template < class K, class V >
class RWTValHashDictionary {

typedef typename hash_dict_t::iterator iterator;
typedef typename hash_dict_t::const_iterator const_iterator;

public:

  RWTValHashDictionary(){}

  RWTValHashDictionary( unsigned (*hashfun)( const K& key) ,size_t n=100):
    fhm(n,G4Hash<K>(hashfun),equal_to<K>()){}

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

typedef typename hash_dict_t::iterator iterator;
typedef typename hash_dict_t::const_iterator const_iterator;
 
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
