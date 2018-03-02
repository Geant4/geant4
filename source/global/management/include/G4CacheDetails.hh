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
//
//    The classes contained in this header files are used by
//    G4Cache to store a TLS instance of the cached object.
//    These classes should not be used by client code, but
//    are used by one of the G4Cache classes.
//
//    G4Cache is a container of the cached value.
//    Not memory efficient, but CPU efficient (constant time access)
//    A different version with a map instead of a vector should be
//    memory efficient and less CPU efficient (log-time access).
//    If really a lot of these objects are used
//    we may want to consider the map version to save some memory
//
//    These are simplified "split-classes" without any
//    copy-from-master logic. Each cached object is associated
//    a unique identified (an integer), that references an instance
//    of the cached value (of template type VALTYPE) in a
//    static TLS data structure
//
//    In case the cache is used for a cached object the object class
//    has to provide a default constructor. Alternatively pointers to
//    objects can be stored in the cache and this limitation is removed
//    but explicit handling of memory (new/delete) of cached object becomes
//    client responsibility.

// History:
//    21 Oct 2013: A. Dotti - First implementation
//
// Todos: - Understand if map based class can be more efficent than
//          vector one
//        - Evaluate use of specialized allocator for TLS "new"
// ------------------------------------------------------------

#ifndef G4CacheDetails_hh
#define G4CacheDetails_hh

#include <vector>
#include "G4Threading.hh"
#include "globals.hh"

#ifdef g4cdebug
  #include <iostream>
  #include <sstream>
  using std::cout;
  using std::endl;
#endif

// A TLS storage for a cache of type VALTYPE
//
template<class VALTYPE> class G4CacheReference
{
  public:

    inline void Initialize( unsigned int id );
    // Initliaze TLS storage

    inline void Destroy( unsigned int id , G4bool last);
    // Cleanup TLS storage for instance id. If last==true
    // destroy and cleanup object container

    inline VALTYPE& GetCache(unsigned int id) const;
    // Returns cached value for instance id

  private:

    typedef std::vector<VALTYPE*> cache_container;
    // Implementation detail: the cached object is stored as a
    // pointer. Object is stored as a pointer to avoid too large
    // std::vector in case of stored objects and allow use of
    // specialized allocators

    static G4ThreadLocal cache_container *cache;
};

// Template specialization for pointers
// Note: Objects are not owned by cache, for this version of the cache
//       the explicit new/delete of the cached object
//
template<class VALTYPE> class G4CacheReference<VALTYPE*>
{
  public:
    inline void Initialize( unsigned int id );
    
    inline void Destroy( unsigned int id , G4bool last);
    
    inline VALTYPE*& GetCache(unsigned int id) const;

  private:
    typedef std::vector<VALTYPE*> cache_container;
    static G4ThreadLocal cache_container *cache;
};

// Template specialization for probably the most used case: double
// Be more efficient avoiding unnecessary "new/delete"
//
template<> class G4CacheReference<G4double>
{
  public:

    inline void Initialize( unsigned int id );
    
    inline void Destroy( unsigned int id , G4bool last);
    
    inline G4double& GetCache(unsigned int id) const;

  private:

    typedef std::vector<G4double> cache_container;
    G4GLOB_DLL static G4ThreadLocal std::vector<G4double> *cache;
};


//================================
// Implementation details follow
//================================

//======= Implementation: G4CacheReference<V>
//===========================================

template<class V>
void G4CacheReference<V>::Initialize( unsigned int id )
{
#ifdef g4cdebug
  if ( cache == 0 )
    cout<<"Generic template"<<endl;
#endif

  // Create cache container
  if ( cache == 0 )
    cache = new cache_container;
  if ( cache->size() <= id )
    cache->resize(id+1,static_cast<V*>(0));
  if ( (*cache)[id] == 0 )
    (*cache)[id]=new V;
}

template<class V>
void G4CacheReference<V>::Destroy( unsigned int id, G4bool last )
{
  if ( cache )
  {
#ifdef g4cdebug
    cout<<"Destroying element"<<id<<" is last?"<<last<<endl;
#endif
    if ( cache->size() < id )
    {
      G4ExceptionDescription msg;
      msg << "Internal fatal error. Invalid G4Cache size (requested id: "
          << id << " but cache has size: "<<cache->size();
      msg << " Possibly client created G4Cache object in a thread and"
          << " tried to delete it from another thread!";
      G4Exception("G4CacheReference<V>::Destroy", "Cache001",
                  FatalException, msg);
      return;
    }
    if ( cache->size() > id && (*cache)[id] )
    {
      delete (*cache)[id];
      (*cache)[id]=0;
    }
    if (last)
    {
      delete cache;
      cache = 0;
    }
  }
}

template<class V>
V& G4CacheReference<V>::GetCache( unsigned int id ) const
{
  return *(cache->operator[](id));
}

template<class V>
G4ThreadLocal typename
G4CacheReference<V>::cache_container * G4CacheReference<V>::cache = 0;

//======= Implementation: G4CacheReference<V*>
//============================================

template<class V>
void G4CacheReference<V*>::Initialize( unsigned int id )
{
#ifdef g4cdebug
  if ( cache == 0 )
    cout<<"Pointer template"<<endl;
#endif
  if ( cache == 0 )
    cache = new cache_container;
  if ( cache->size() <= id )
    cache->resize(id+1,static_cast<V*>(0));
}

template<class V>
inline void G4CacheReference<V*>::Destroy( unsigned int id , G4bool last )
{
  if ( cache )
  {
#ifdef g4cdebug
    cout << "Destroying element" << id << " is last?" << last
         << "-Pointer template specialization-" << endl;
#endif
    if ( cache->size() < id )
    {
      G4ExceptionDescription msg;
      msg << "Internal fatal error. Invalid G4Cache size (requested id: "
          << id << " but cache has size: " << cache->size();
      msg << " Possibly client created G4Cache object in a thread and"
          << " tried to delete it from another thread!";
      G4Exception("G4CacheReference<V*>::Destroy", "Cache001",
                  FatalException, msg);
      return;
    }
    if ( cache->size() > id && (*cache)[id] )
    {
      // Ownership is for client
      // delete (*cache)[id];
      (*cache)[id]=0;
    }
    if (last )
    {
      delete cache;
      cache = 0;
    }
  }
}

template<class V>
V*& G4CacheReference<V*>::GetCache(unsigned int id) const
{
  return (cache->operator[](id));
}

template<class V>
G4ThreadLocal typename
G4CacheReference<V*>::cache_container * G4CacheReference<V*>::cache = 0;

//======= Implementation: G4CacheReference<double>
//============================================

void G4CacheReference<G4double>::Initialize( unsigned int id )
{
#ifdef g4cdebug
  cout<<"Specialized template for G4double"<<endl;
#endif
  if ( cache == 0 )
    cache = new cache_container;
  if ( cache->size() <= id )
    cache->resize(id+1,static_cast<G4double>(0));
}

#ifdef g4cdebug
void G4CacheReference<G4double>::Destroy( unsigned int id , G4bool last)
#else
void G4CacheReference<G4double>::Destroy( unsigned int /*id*/ , G4bool last)
#endif
{
  if ( cache && last )
  {
#ifdef g4cdebug
    cout << "Destroying element" << id << " is last?" << last
         << "-Pointer template specialization-" << endl;
#endif
    delete cache;
    cache = 0;
  }
}

G4double& G4CacheReference<G4double>::GetCache(unsigned int id) const
{
  return cache->operator[](id);
}

#endif
