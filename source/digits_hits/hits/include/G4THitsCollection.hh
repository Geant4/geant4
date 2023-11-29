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
//
//

#ifndef G4THitsCollection_h
#define G4THitsCollection_h 1

#include "G4Allocator.hh"
#include "G4VHitsCollection.hh"
#include "globals.hh"

#include <vector>

// class description:
//
//  This is a template class of hits collection and parametrized by
// The concrete class of G4VHit. This is a uniform collection for
// a particular concrete hit class objects.
//  An intermediate layer class G4HitsCollection appeared in this
// header file is used just for G4Allocator, because G4Allocator
// cannot be instansiated with a template class. Thus G4HitsCollection
// class MUST NOT be directly used by the user.

class G4HitsCollection : public G4VHitsCollection
{
 public:
  using G4VHitsCollection::G4VHitsCollection;
  ~G4HitsCollection() override = default;
  G4bool operator==(const G4HitsCollection& right) const;

 protected:
  void* theCollection = nullptr;
};

#if defined G4DIGI_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4HitsCollection>*& anHCAllocator_G4MT_TLS_();
#else
extern G4DLLIMPORT G4Allocator<G4HitsCollection>*& anHCAllocator_G4MT_TLS_();
#endif

template <class T>
class G4THitsCollection : public G4HitsCollection
{
 public:
  G4THitsCollection();
  G4THitsCollection(G4String detName, G4String colNam);
  ~G4THitsCollection() override;

  G4bool operator==(const G4THitsCollection<T>& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* anHC);

  // Invoke Draw() method on all hit objects in collection
  void DrawAllHits() override;

  // Invoke Print() method on all hit objects in collection
  void PrintAllHits() override;

  //  Returns a pointer to the concrete hit object
  // Returns pointer to the concrete hit object at given index
  // Not bounds checked
  inline T* operator[](size_t i) const { return (*((std::vector<T*>*)theCollection))[i]; }

  // Return pointer to hit collection
  inline std::vector<T*>* GetVector() const { return (std::vector<T*>*)theCollection; }

  // Insert a hit object in the collection, taking onwenership
  // Returns the total number of hit objects stored after insertion
  inline size_t insert(T* aHit)
  {
    auto theHitsCollection = (std::vector<T*>*)theCollection;
    theHitsCollection->push_back(aHit);
    return theHitsCollection->size();
  }

  // Returns the number of hit objects stored in this collection.
  inline size_t entries() const
  {
    auto theHitsCollection = (std::vector<T*>*)theCollection;
    return theHitsCollection->size();
  }

  G4VHit* GetHit(size_t i) const override { return (*((std::vector<T*>*)theCollection))[i]; }

  size_t GetSize() const override { return ((std::vector<T*>*)theCollection)->size(); }
};

template <class T>
inline void* G4THitsCollection<T>::operator new(size_t)
{
  if (anHCAllocator_G4MT_TLS_() == nullptr) {
    anHCAllocator_G4MT_TLS_() = new G4Allocator<G4HitsCollection>;
  }
  return (void*)anHCAllocator_G4MT_TLS_()->MallocSingle();
}

template <class T>
inline void G4THitsCollection<T>::operator delete(void* anHC)
{
  anHCAllocator_G4MT_TLS_()->FreeSingle((G4HitsCollection*)anHC);
}

template <class T>
G4THitsCollection<T>::G4THitsCollection()
{
  auto theHitsCollection = new std::vector<T*>;
  theCollection = (void*)theHitsCollection;
}

template <class T>
G4THitsCollection<T>::G4THitsCollection(G4String detName, G4String colNam)
  : G4HitsCollection(detName, colNam)
{
  auto theHitsCollection = new std::vector<T*>;
  theCollection = (void*)theHitsCollection;
}

template <class T>
G4THitsCollection<T>::~G4THitsCollection()
{
  auto theHitsCollection = (std::vector<T*>*)theCollection;
  for (const auto* hit : *theHitsCollection) {
    delete hit;
  }
  theHitsCollection->clear();
  delete theHitsCollection;
}

template <class T>
G4bool G4THitsCollection<T>::operator==(const G4THitsCollection<T>& right) const
{
  return (collectionName == right.collectionName);
}

template <class T>
void G4THitsCollection<T>::DrawAllHits()
{
  auto theHitsCollection = (std::vector<T*>*)theCollection;
  for (auto* hit : *theHitsCollection) {
    hit->Draw();
  }
}

template <class T>
void G4THitsCollection<T>::PrintAllHits()
{
  auto theHitsCollection = (std::vector<T*>*)theCollection;
  for (auto* hit : *theHitsCollection) {
    hit->Print();
  }
}

#endif
