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
// G4CacheDetails
//
// Class description:
//
//    The classes contained in this header files are used by
//    G4Cache to store a TLS instance of the cached object.
//    These classes should not be used by client code, but
//    are used by one of the G4Cache classes.
//
//    G4Cache is a container of the cached value.
//    Not memory efficient, but CPU efficient (constant time access).
//    A different version with a map instead of a vector should be
//    memory efficient and less CPU efficient (log-time access).
//    If really a lot of these objects are used
//    we may want to consider the map version to save some memory.
//
//    These are simplified "split-classes" without any
//    copy-from-master logic. Each cached object is associated
//    a unique identified (an integer), that references an instance
//    of the cached value (of template type VALTYPE) in a
//    static TLS data structure.
//
//    In case the cache is used for a cached object the object class
//    has to provide a default constructor. Alternatively pointers to
//    objects can be stored in the cache and this limitation is removed
//    but explicit handling of memory (new/delete) of cached object becomes
//    client responsibility.
//
// Todos: - Understand if map based class can be more efficent than
//          vector one.
//        - Evaluate use of specialized allocator for TLS "new".

// Author: A.Dotti, 21 October 2013 - First implementation
// --------------------------------------------------------------------
#ifndef G4CacheDetails_hh
#define G4CacheDetails_hh

#include "G4Threading.hh"
#include "globals.hh"
#include <vector>

// A TLS storage for a cache of type VALTYPE
//
template <class VALTYPE>
class G4CacheReference
{
 public:
  inline void Initialize(unsigned int id);
  // Initliaze TLS storage

  inline void Destroy(unsigned int id, G4bool last);
  // Cleanup TLS storage for instance id. If last==true
  // destroy and cleanup object container

  inline VALTYPE& GetCache(unsigned int id) const;
  // Returns cached value for instance id

 private:
  using cache_container = std::vector<VALTYPE*>;
  // Implementation detail: the cached object is stored as a
  // pointer. Object is stored as a pointer to avoid too large
  // std::vector in case of stored objects and allow use of
  // specialized allocators

  static cache_container*& cache();
};

// Template specialization for pointers
// Note: Objects are not owned by cache, for this version of the cache
//       the explicit new/delete of the cached object
//
template <class VALTYPE>
class G4CacheReference<VALTYPE*>
{
 public:
  inline void Initialize(unsigned int id);

  inline void Destroy(unsigned int id, G4bool last);

  inline VALTYPE*& GetCache(unsigned int id) const;

 private:
  using cache_container = std::vector<VALTYPE*>;
  static cache_container*& cache();
};

// Template specialization for probably the most used case: double
// Be more efficient avoiding unnecessary "new/delete"
//
template <>
class G4CacheReference<G4double>
{
 public:
  inline void Initialize(unsigned int id);

  inline void Destroy(unsigned int id, G4bool last);

  inline G4double& GetCache(unsigned int id) const;

 private:
  using cache_container = std::vector<G4double>;
  static G4GLOB_DLL cache_container*& cache();
};

//======= Implementation: G4CacheReference<V>
//===========================================

template <class V>
void G4CacheReference<V>::Initialize(unsigned int id)
{
  // Create cache container
  if(cache() == nullptr)
  {
#ifdef g4cdebug
    std::cout << "Generic template container..." << std::endl;
#endif
    cache() = new cache_container;
  }
  if(cache()->size() <= id)
  {
    cache()->resize(id + 1, static_cast<V*>(0));
  }
  if((*cache())[id] == 0)
  {
    (*cache())[id] = new V;
  }
}

template <class V>
void G4CacheReference<V>::Destroy(unsigned int id, G4bool last)
{
  if(cache() != nullptr)
  {
#ifdef g4cdebug
    std::cout << "V: Destroying element " << id << " is last? " << last
              << std::endl;
#endif
    if(cache()->size() < id)
    {
      G4ExceptionDescription msg;
      msg << "Internal fatal error. Invalid G4Cache size (requested id: " << id
          << " but cache has size: " << cache()->size();
      msg << " Possibly client created G4Cache object in a thread and"
          << " tried to delete it from another thread!";
      G4Exception("G4CacheReference<V>::Destroy", "Cache001", FatalException,
                  msg);
      return;
    }
    if(cache()->size() > id && (*cache())[id] != nullptr)
    {
#ifdef g4cdebug
      std::cout << "V: Destroying element " << id
                << " size: " << cache()->size() << std::endl;
#endif
      delete(*cache())[id];
      (*cache())[id] = nullptr;
    }
    if(last)
    {
#ifdef g4cdebug
      std::cout << "V: Destroying LAST element!" << std::endl;
#endif
      delete cache();
      cache() = nullptr;
    }
  }
}

template <class V>
V& G4CacheReference<V>::GetCache(unsigned int id) const
{
  return *(cache()->operator[](id));
}

template <class V>
typename G4CacheReference<V>::cache_container*& G4CacheReference<V>::cache()
{
  G4ThreadLocalStatic cache_container* _instance = nullptr;
  return _instance;
}

//======= Implementation: G4CacheReference<V*>
//============================================

template <class V>
void G4CacheReference<V*>::Initialize(unsigned int id)
{
  if(cache() == nullptr)
  {
#ifdef g4cdebug
    std::cout << "Pointer template container..." << std::endl;
#endif
    cache() = new cache_container;
  }
  if(cache()->size() <= id)
  {
    cache()->resize(id + 1, static_cast<V*>(nullptr));
  }
}

template <class V>
inline void G4CacheReference<V*>::Destroy(unsigned int id, G4bool last)
{
  if(cache() != nullptr)
  {
#ifdef g4cdebug
    std::cout << "V*: Destroying element " << id << " is last? " << last
              << std::endl;
#endif
    if(cache()->size() < id)
    {
      G4ExceptionDescription msg;
      msg << "Internal fatal error. Invalid G4Cache size (requested id: " << id
          << " but cache has size: " << cache()->size();
      msg << " Possibly client created G4Cache object in a thread and"
          << " tried to delete it from another thread!";
      G4Exception("G4CacheReference<V*>::Destroy", "Cache001", FatalException,
                  msg);
      return;
    }
    if(cache()->size() > id && (*cache())[id] != nullptr)
    {
      // Ownership is for client
      // delete (*cache())[id];
#ifdef g4cdebug
      std::cout << "V*: Resetting element " << id
                << " size: " << cache()->size() << std::endl;
#endif
      (*cache())[id] = nullptr;
    }
    if(last)
    {
#ifdef g4cdebug
      std::cout << "V*: Deleting LAST element!" << std::endl;
#endif
      delete cache();
      cache() = nullptr;
    }
  }
}

template <class V>
V*& G4CacheReference<V*>::GetCache(unsigned int id) const
{
  return (cache()->operator[](id));
}

template <class V>
typename G4CacheReference<V*>::cache_container*& G4CacheReference<V*>::cache()
{
  G4ThreadLocalStatic cache_container* _instance = nullptr;
  return _instance;
}

//======= Implementation: G4CacheReference<double>
//============================================

void G4CacheReference<G4double>::Initialize(unsigned int id)
{
  if(cache() == nullptr)
  {
#ifdef g4cdebug
    std::cout << "Specialized template for G4double container..." << std::endl;
#endif
    cache() = new cache_container;
  }
  if(cache()->size() <= id)
  {
    cache()->resize(id + 1, static_cast<G4double>(0));
  }
}

void G4CacheReference<G4double>::Destroy(unsigned int /*id*/, G4bool last)
{
  if(cache() != nullptr && last)
  {
#ifdef g4cdebug
    std::cout << "DB: Destroying LAST element! Is it last? " << last
              << std::endl;
#endif
    delete cache();
    cache() = nullptr;
  }
}

G4double& G4CacheReference<G4double>::GetCache(unsigned int id) const
{
  return cache()->operator[](id);
}

#endif
