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
/// \file parallel/ThreadsafeScorers/include/G4TAtomicHitsMap.hh
/// \brief Definition of the G4TAtomicHitsMap class
//
//
//
//
/// This is an implementation of G4THitsMap<T> where the underlying
///     type is G4atomic<T>, not just T. A static assert is provided to
///     ensure that T is fundamental. This class should be used in lieu
///     of G4THitsMap<T> when memory is a concern. Atomics are
///     thread-safe and *generally* faster that mutexes (as long as the
///     STL implementation is lock-free) but the synchronization does
///     not come without a cost. If performance is the primary concern,
///     use G4THitsMap<T> in thread-local instances.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




#ifndef G4TAtomicHitsMap_h
#define G4TAtomicHitsMap_h 1

#include "G4THitsCollection.hh"
#include "G4THitsMap.hh"
#include "globals.hh"
#include "G4atomic.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <map>
#include <type_traits>

// class description:
//
//  This is a template class of hits map and parametrized by
// The concrete class of G4VHit. This is a uniform collection for
// a particular concrete hit class objects.
//  An intermediate layer class G4HitsMap appeared in this
// header file is used just for G4Allocator, because G4Allocator
// cannot be instansiated with a template class. Thus G4HitsMap
// class MUST NOT be directly used by the user.

template <typename T>
class G4TAtomicHitsMap : public G4VHitsCollection
{
protected:
  static_assert(std::is_fundamental<T>::value,
                "G4TAtomicHitsMap must use fundamental type");

public:
  typedef G4atomic<T>                             value_type;
  typedef value_type*                             mapped_type;
  typedef typename std::map<G4int, mapped_type>   container_type;
  typedef typename container_type::iterator       iterator;
  typedef typename container_type::const_iterator const_iterator;

public:
  G4TAtomicHitsMap();

public: // with description
  G4TAtomicHitsMap(G4String detName, G4String colNam);
  // constructor.

public:
  virtual ~G4TAtomicHitsMap();
  G4bool operator==(const G4TAtomicHitsMap<T> &right) const;
  G4TAtomicHitsMap<T> & operator+=(const G4TAtomicHitsMap<T> &right) const;
  G4TAtomicHitsMap<T> & operator+=(const G4THitsMap<T> &right) const;

public: // with description
  virtual void DrawAllHits();
  virtual void PrintAllHits();
  //  These two methods invokes Draw() and Print() methods of all of
  // hit objects stored in this map, respectively.

public: // with description
  inline value_type* operator[](G4int key) const;

  //  Returns a pointer to a concrete hit object.
  inline container_type* GetMap() const
  { return theCollection; }
  //  Returns a collection map.
  inline G4int add(const G4int & key, value_type*& aHit) const;
  inline G4int add(const G4int & key, T& aHit) const;
  //  Insert a hit object. Total number of hit objects stored in this
  // map is returned.
  inline G4int set(const G4int & key, value_type*& aHit) const;
  inline G4int set(const G4int & key, T& aHit) const;
  //  Overwrite a hit object. Total number of hit objects stored in this
  // map is returned.
  inline G4int entries() const
  {
    return theCollection->size();
  }
  //  Returns the number of hit objects stored in this map
  inline void clear();

public:
  virtual G4VHit* GetHit(size_t) const {return 0;}
  virtual size_t GetSize() const
  {
    return theCollection->size();
  }

  virtual size_t size() const { return theCollection->size(); }

public:
  iterator begin() { return theCollection->begin(); }
  iterator end() { return theCollection->end(); }

  const_iterator begin() const { return theCollection->begin(); }
  const_iterator end() const { return theCollection->end(); }

  const_iterator cbegin() const { return theCollection->cbegin(); }
  const_iterator cend() const { return theCollection->cend(); }

  iterator find(G4int p) { return theCollection->find(p); }
  const_iterator find(G4int p) const { return theCollection->find(p); }

private:
  container_type*   theCollection;
  mutable G4Mutex   fMutex;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
G4TAtomicHitsMap<T>::G4TAtomicHitsMap()
  : theCollection(new container_type)
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
G4TAtomicHitsMap<T>::G4TAtomicHitsMap(G4String detName,
                                      G4String colNam)
  : G4VHitsCollection(detName,colNam),
    theCollection(new container_type)
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
G4TAtomicHitsMap<T>::~G4TAtomicHitsMap()
{
  for(auto itr = theCollection->begin(); itr != theCollection->end(); itr++)
      delete itr->second;

  delete theCollection;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
G4bool G4TAtomicHitsMap<T>::operator==(const G4TAtomicHitsMap<T> &right) const
{
  return (collectionName == right.collectionName);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
G4TAtomicHitsMap<T>&
G4TAtomicHitsMap<T>::operator+=(const G4TAtomicHitsMap<T>& rhs) const
{
  for(auto itr = rhs.GetMap()->begin(); itr != rhs.GetMap()->end(); itr++)
    add(itr->first, *(itr->second));

  return (G4TAtomicHitsMap<T>&)(*this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
G4TAtomicHitsMap<T>&
G4TAtomicHitsMap<T>::operator+=(const G4THitsMap<T>& rhs) const
{
  for(auto itr = rhs.GetMap()->begin(); itr != rhs.GetMap()->end(); itr++)
    add(itr->first, *(itr->second));

  return (G4TAtomicHitsMap<T>&)(*this);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
inline typename G4TAtomicHitsMap<T>::value_type*
G4TAtomicHitsMap<T>::operator[](G4int key) const
{
  if(theCollection->find(key) != theCollection->end())
    return theCollection->find(key)->second;
  else
  {
    G4AutoLock l(&fMutex);
    if(theCollection->find(key) == theCollection->end())
    {
      value_type* ptr = new value_type;
      (*theCollection)[key] = ptr;
      return ptr;
    } else
      return theCollection->find(key)->second;
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
inline G4int
G4TAtomicHitsMap<T>::add(const G4int& key, value_type*& aHit) const
{
  if(theCollection->find(key) != theCollection->end())
    *(*theCollection)[key] += *aHit;
  else
  {
    G4AutoLock l(&fMutex);
    (*theCollection)[key] = aHit;
  }
  G4AutoLock l(&fMutex);
  return theCollection->size();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
inline G4int
G4TAtomicHitsMap<T>::add(const G4int& key, T& aHit) const
{

  if(theCollection->find(key) != theCollection->end())
    *(*theCollection)[key] += aHit;
  else
  {
    value_type* hit = new value_type;
    *hit = aHit;
    G4AutoLock l(&fMutex);
    (*theCollection)[key] = hit;
  }
  G4AutoLock l(&fMutex);
  return theCollection->size();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
inline G4int
G4TAtomicHitsMap<T>::set(const G4int& key, value_type*& aHit) const
{
  if(theCollection->find(key) != theCollection->end())
      delete (*theCollection)[key]->second;

  (*theCollection)[key] = aHit;
  G4AutoLock l(&fMutex);
  return theCollection->size();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
inline G4int
G4TAtomicHitsMap<T>::set(const G4int& key, T& aHit) const
{
    if(theCollection->find(key) != theCollection->end())
        *(*theCollection)[key] = aHit;
    else
    {
      value_type* hit = new value_type;
      *hit = aHit;
      (*theCollection)[key] = hit;
    }
    G4AutoLock l(&fMutex);
    return theCollection->size();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
void G4TAtomicHitsMap<T>::DrawAllHits()
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
void G4TAtomicHitsMap<T>::PrintAllHits()
{
  G4cout << "G4TAtomicHitsMap " << SDname << " / " << collectionName << " --- "
         << entries() << " entries" << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <typename T>
void G4TAtomicHitsMap<T>::clear()
{
  G4AutoLock l(&fMutex);

  for(auto itr = theCollection->begin(); itr != theCollection->end(); itr++)
    delete itr->second;

  theCollection->clear();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
