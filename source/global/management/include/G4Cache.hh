//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
// ---------------------------------------------------------------
// GEANT 4 class header file
//
// Class Description:
//      Helper classes for Geant4 Multi-Threaded.
//      The classes defined in this header file provide a thread-private
//      cache to store, in a class instance shared among threads,
//      a thread-local variable V.
//      These are templated classes on the to-be-stored object.
//
// Example:
//      Let's assume an instance myObject of class G4Shared is sharead between
//      threads. Still a data member of this class needs to be thread-private.
//      A typical example of this being a "cache" for a local calculation.
//      The helper here defined can be used to guarantee thread-safe operations
//      on the thread-private object.
//      Example:
//      class G4Shared {
//          G4double sharedData;
//          G4Cache<G4double> threadPrivate;
//          void foo() {
//              G4double priv = threadPrivate.Get();
//              if ( priv < 10 ) priv += sharedData;
//              threadPrivate.Put( priv );
//          }
//      };
//
//      Two variants of the base G4Cache exists. The first one being
//      G4VectorCache similar to std::vector
//      Example:
//          G4VectorCache<G4double> aVect;
//          aVect.Push_back( 3.2 );
//          aVect.Push_back( 4.1 );
//          cout<<aVect[0]<<endl;
//      The second one being:
//      G4MapCache    similar to std::map
//      Example:
//          G4MapCache<G4int,G4double> aMap;
//          aMap[320]=1.234;
//
//      See classes definition for details.
//      See testG4Cache unit test for details on usage
//
// History:
//  21 October 2013: A. Dotti - First implementation

#ifndef G4CACHE_HH
#define G4CACHE_HH

//Debug this code
//#define g4cdebug 1

//Thread Local storage details are in this header file
#include "G4CacheDetails.hh"

// A templated cache to store a thread-private data of type VALTYPE.
template<class VALTYPE>
class G4Cache {
public:
    typedef VALTYPE value_type;
    // The stored type
  
    G4Cache();
    // Default constructor
    
    G4Cache(const value_type& v);
    // Construct cache object with initial value

    virtual ~G4Cache();
    // Default destructor
    
    inline value_type& Get() const;
    // Gets reference to cached value of this threads
  
    inline void Put( const value_type& val ) const;
    // Sets this thread cached value to val

    inline value_type Pop();
    // Gets copy of cached value
    
    G4Cache(const G4Cache& rhs);
    G4Cache& operator=(const G4Cache& rhs);

protected:
  const int& GetId() const { return id; }
private:
  int id;
  mutable G4CacheReference<value_type> theCache;
  static G4Mutex gMutex;
  static unsigned int instancesctr;
  static unsigned int dstrctr;

  inline value_type& GetCache() const {
    theCache.Initialize(id);
    return theCache.GetCache(id);
  }

};


// A vector version of the cache. Implements vector interface.
// Can be used directly as a std::vector would be used.
template<class VALTYPE>
class G4VectorCache : public G4Cache< std::vector<VALTYPE> > {
public:
    //Some useful defintions
    typedef VALTYPE value_type;
    typedef typename std::vector<value_type> vector_type;
    typedef typename vector_type::size_type size_type;
    typedef typename vector_type::iterator iterator;
    typedef typename vector_type::const_iterator const_iterator;

    G4VectorCache();
    // Default constructor
    
    G4VectorCache( G4int nElems );
    // Creates a vector cache of nElems elements
    
    G4VectorCache( G4int nElems , value_type* vals );
    // Creates a vector cache with elements from an array
    
    virtual ~G4VectorCache();
    // Default destructor

    // Interface with funxtionalities of similar name of std::vector
    inline void Push_back( const value_type& val );
    inline value_type Pop_back();
    inline value_type& operator[](const G4int& idx);
    inline iterator Begin();
    inline iterator End();
    inline void Clear();
    inline size_type Size() { return G4Cache<vector_type>::Get().size(); } //Needs to be here for a VC9 compilation problem
};


// a Map version of the cache. Implemetns std::map interface.
// Can be used directly as a std::map would be used.
// KEYTYPE being the key type and VALTYPE the value type.
#include <map>
template<class KEYTYPE, class VALTYPE>
class G4MapCache : public G4Cache<std::map<KEYTYPE,VALTYPE> > {
public:
    //Some useful definitions
    typedef KEYTYPE key_type;
    typedef VALTYPE value_type;
    typedef typename std::map<key_type,value_type> map_type;
    typedef typename map_type::size_type size_type;
    typedef typename map_type::iterator iterator;
    typedef typename map_type::const_iterator const_iterator;

    virtual ~G4MapCache();
    // Default destructor

    inline G4bool Has(const key_type& k );
    // Returns true if map contains element corresponding to key k
    
    // Interface with functionalities of similar name of std::map
    inline std::pair<iterator,G4bool> Insert( const key_type& k , const value_type& v );
    inline iterator Begin();
    inline iterator End();
    inline iterator Find(const key_type& k );
    inline value_type& Get(const key_type& k );
    inline size_type Erase(const key_type& k );
    inline value_type& operator[](const key_type& k);
    inline size_type Size() { return G4Cache<map_type>::Get().size(); } //Needs to be here for a VC9 compilation problem
};



//=============================================================
// Implementation details follow
//=============================================================



#ifdef g4cdebug
#include <iostream>
#include <sstream>
using std::cout;
using std::endl;
#endif

#include "G4AutoLock.hh"


//========= Implementation: G4Cache<V>

template<class V>
G4Cache<V>::G4Cache()
{
    G4AutoLock l(&gMutex);
    id = instancesctr++;
#ifdef g4cdebug
    cout<<"G4Cache id: "<<id<<endl;
#endif
}

template<class V>
G4Cache<V>::G4Cache(const G4Cache<V>& rhs)
{
	//Copy is special, we need to copy the content
	//of the cache, not the cache object
	if ( this == &rhs ) return;
	G4AutoLock l(&gMutex);
	id = instancesctr++;
	//Force copy of cached data
	V aCopy = rhs.GetCache();
	Put( aCopy );
#ifdef g4cdebug
	cout<<"Copy constructor with id: "<<id<<endl;
#endif
}

template<class V>
G4Cache<V>& G4Cache<V>::operator=(const G4Cache<V>& rhs)
{
	if (this == &rhs) return *this;
	//Force copy of cached data
	V aCopy = rhs.GetCache();
	Put(aCopy);
#ifdef g4cdebug
	cout<<"Assignement operator with id: "<<id<<endl;
#endif
	return *this;
}

template<class V>
G4Cache<V>::G4Cache(const V& v)
{
    G4AutoLock l(&gMutex);
    id = instancesctr++;
    Put(v);
#ifdef g4cdebug
    cout<<"G4Cache id: "<<id<<" "<<endl;
#endif
}

template<class V>
G4Cache<V>::~G4Cache()
{ //Move base calss
#ifdef g4cdebug
    cout<<"~G4Cache id: "<<id<<" "<<endl;
#endif
    G4AutoLock l(&gMutex);
    ++dstrctr;
    G4bool last = ( dstrctr == instancesctr );
    theCache.Destroy(id,last);
    if (last) {
        instancesctr = 0;
        dstrctr = 0;
    }
}

template<class V>
V& G4Cache<V>::Get() const
{ return GetCache(); }

template<class V>
void G4Cache<V>::Put( const V& val ) const
{ GetCache() = val; }

//Should here remove from cache element?
template<class V>
V G4Cache<V>::Pop()
{ return GetCache(); }

template<class V>
unsigned int G4Cache<V>::instancesctr = 0;

template<class V>
unsigned int G4Cache<V>::dstrctr = 0;

template<class V>
G4Mutex G4Cache<V>::gMutex = G4MUTEX_INITIALIZER;

//========== Implementation: G4VectorCache<V>
template<class V>
G4VectorCache<V>::G4VectorCache()
{ }

template<class V>
G4VectorCache<V>::~G4VectorCache() {
#ifdef g4cdebug
    cout<<"~G4VectorCache "<<G4Cache<G4VectorCache<V>::vector_type>::GetId()<<" with size: "<<Size()<<"->";
    for ( size_type i = 0 ; i < Size() ; ++i )
        cout<<operator[](i)<<",";
    cout<<"<-"<<endl;
#endif
}

template<class V>
G4VectorCache<V>::G4VectorCache(G4int nElems ) {
    vector_type& cc = G4Cache<vector_type>::Get();
    cc.resize(nElems);
}

template<class V>
G4VectorCache<V>::G4VectorCache(G4int nElems , V* vals ) {
    vector_type& cc = G4Cache<vector_type>::Get();
    cc.resize(nElems);
    for ( G4int idx = 0 ; idx < nElems ; ++idx )
        cc[idx]=vals[idx];
}

template<class V>
void G4VectorCache<V>::Push_back( const V& val )
{
    G4Cache<vector_type>::Get().push_back( val );
}

template<class V>
V G4VectorCache<V>::Pop_back()
{
    vector_type& cc = G4Cache<vector_type>::Get();
    value_type val = cc[cc.size()-1];
    cc.pop_back();
    return val;
}

template<class V>
V& G4VectorCache<V>::operator[](const G4int& idx)
{
    vector_type& cc = G4Cache<vector_type>::Get();
    return cc[idx];
}

template<class V>
typename G4VectorCache<V>::iterator G4VectorCache<V>::Begin()
{
    return G4Cache<vector_type>::Get().begin();
}

template<class V>
typename G4VectorCache<V>::iterator G4VectorCache<V>::End()
{
    return G4Cache<vector_type>::Get().end();
}

template<class V>
void G4VectorCache<V>::Clear()
{
    G4Cache<vector_type>::Get().clear();
}

//template<class V>
//typename G4VectorCache<V>::size_type G4VectorCache<V>::Size()
//{
//    return G4Cache<vector_type>::Get().size();
//}

//======== Implementation: G4MapType<K,V>
template<class K, class V>
G4MapCache<K,V>::~G4MapCache()
{
#ifdef g4cdebug
    cout<<"~G4MacCache "<<G4Cache<map_type>::GetId()<<" with size: "<<Size()<<"->";
    for ( iterator it = Begin() ; it != End() ; ++it )
        cout<<it->first<<":"<<it->second<<",";
    cout<<"<-"<<endl;
#endif
}

template<class K, class V>
std::pair<typename G4MapCache<K,V>::iterator,G4bool> G4MapCache<K,V>::Insert(
                                                                            const K& k,
                                                                            const V& v
                                                                             )
{
    return G4Cache<map_type>::Get().insert( std::pair<key_type,value_type>(k,v) );
}

//template<class K, class V>
//typename G4MapCache<K,V>::size_type G4MapCache<K,V>::Size()
//{
//    return G4Cache<map_type>::Get().size();
//}

template<class K, class V>
typename G4MapCache<K,V>::iterator G4MapCache<K,V>::Begin()
{
    return G4Cache<map_type>::Get().begin();
}
template<class K, class V>
typename G4MapCache<K,V>::iterator G4MapCache<K,V>::End()
{
    return G4Cache<map_type>::Get().end();
}

template<class K, class V>
typename G4MapCache<K,V>::iterator G4MapCache<K,V>::Find(const K& k )
{
    return G4Cache<map_type>::Get().find(k);
}

template<class K, class V>
G4bool G4MapCache<K,V>::Has(const K& k )
{
    return ( Find(k) != End() );
}

template<class K, class V>
V& G4MapCache<K,V>::Get(const K& k )
{
    return Find(k)->second;
}

template<class K, class V>
typename G4MapCache<K,V>::size_type G4MapCache<K,V>::Erase(const K& k )
{
    return G4Cache<map_type>::Get().erase(k);
}

template<class K, class V>
V& G4MapCache<K,V>::operator[](const K& k)
{
    return (G4Cache<map_type>::Get())[k];
}

#endif
