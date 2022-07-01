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

#ifndef G4TDigiCollection_h
#define G4TDigiCollection_h 1

#include "G4VDigiCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
#include <vector>

// class description:
//
//  This is a template class of digi collection and parametrized by
// The concrete class of G4VDigi. This is a uniform collection for
// a particular concrete digi class objects.
//  An intermediate layer class G4DigiCollection appeared in this
// header file is used just for G4Allocator, because G4Allocator
// cannot be instansiated with a template class. Thus G4DigiCollection
// class MUST NOT be directly used by the user.

class G4DigiCollection : public G4VDigiCollection
{
 public:
  G4DigiCollection();
  G4DigiCollection(G4String detName, G4String colNam);
  virtual ~G4DigiCollection();
  G4bool operator==(const G4DigiCollection& right) const;

 protected:
  void* theCollection;
};

#if defined G4DIGI_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4DigiCollection>*& aDCAllocator_G4MT_TLS_();
#else
extern G4DLLIMPORT G4Allocator<G4DigiCollection>*& aDCAllocator_G4MT_TLS_();
#endif

template <class T>
class G4TDigiCollection : public G4DigiCollection
{
 public:
  G4TDigiCollection();

 public:  // with description
  G4TDigiCollection(G4String detName, G4String colNam);
  // Constructor.
 public:
  virtual ~G4TDigiCollection();
  G4bool operator==(const G4TDigiCollection& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aDC);

 public:  // with description
  virtual void DrawAllDigi();
  virtual void PrintAllDigi();
  //  These two methods invokes Draw() and Print() methods of all of
  // digit objects stored in this collection, respectively.

 public:  // with description
  inline T* operator[](size_t i) const
  {
    return (*((std::vector<T*>*) theCollection))[i];
  }
  //  Returns a pointer to a concrete digi object.
  inline std::vector<T*>* GetVector() const
  {
    return (std::vector<T*>*) theCollection;
  }
  //  Returns a collection vector.
  inline size_t insert(T* aHit)
  {
    std::vector<T*>* theDigiCollection = (std::vector<T*>*) theCollection;
    theDigiCollection->push_back(aHit);
    return theDigiCollection->size();
  }
  //  Insert a digi object. Total number of digi objects stored in this
  // collection is returned.
  inline size_t entries() const
  {
    std::vector<T*>* theDigiCollection = (std::vector<T*>*) theCollection;
    return theDigiCollection->size();
  }
  //  Returns the number of digi objcets stored in this collection.

 public:
  virtual G4VDigi* GetDigi(size_t i) const
  {
    return (*((std::vector<T*>*) theCollection))[i];
  }
  virtual size_t GetSize() const
  {
    return ((std::vector<T*>*) theCollection)->size();
  }
};

template <class T>
inline void* G4TDigiCollection<T>::operator new(size_t)
{
  if(aDCAllocator_G4MT_TLS_() == nullptr)
  {
    aDCAllocator_G4MT_TLS_() = new G4Allocator<G4DigiCollection>;
  }
  return (void*) aDCAllocator_G4MT_TLS_()->MallocSingle();
}

template <class T>
inline void G4TDigiCollection<T>::operator delete(void* aDC)
{
  aDCAllocator_G4MT_TLS_()->FreeSingle((G4DigiCollection*) aDC);
}

template <class T>
G4TDigiCollection<T>::G4TDigiCollection()
{
  std::vector<T*>* theDigiCollection = new std::vector<T*>;
  theCollection                      = (void*) theDigiCollection;
}

template <class T>
G4TDigiCollection<T>::G4TDigiCollection(G4String detName, G4String colNam)
  : G4DigiCollection(detName, colNam)
{
  std::vector<T*>* theDigiCollection = new std::vector<T*>;
  theCollection                      = (void*) theDigiCollection;
}

template <class T>
G4TDigiCollection<T>::~G4TDigiCollection()
{
  std::vector<T*>* theDigiCollection = (std::vector<T*>*) theCollection;
  
  for(size_t i = 0; i < theDigiCollection->size(); i++)
  {
    delete(*theDigiCollection)[i];
  }
  theDigiCollection->clear();
  delete theDigiCollection;
}

template <class T>
G4bool G4TDigiCollection<T>::operator==(const G4TDigiCollection<T>& right) const
{
  return (collectionName == right.collectionName);
}

template <class T>
void G4TDigiCollection<T>::DrawAllDigi()
{
  std::vector<T*>* theDigiCollection = (std::vector<T*>*) theCollection;
  size_t n                           = theDigiCollection->size();
  for(size_t i = 0; i < n; i++)
  {
    (*theDigiCollection)[i]->Draw();
  }
}

template <class T>
void G4TDigiCollection<T>::PrintAllDigi()
{
  std::vector<T*>* theDigiCollection = (std::vector<T*>*) theCollection;
  size_t n                           = theDigiCollection->size();
  for(size_t i = 0; i < n; i++)
  {
    (*theDigiCollection)[i]->Print();
  }
}

#endif
