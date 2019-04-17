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
/// \file parallel/ThreadsafeScorers/include/G4TAtomicHitsCollection.hh
/// \brief Definition of the G4TAtomicHitsCollection class
//
//
//
//
/// This is an implementation of G4THitsCollection<T> where the underlying
///     type is G4atomic<T>, not just T. A static assert is provided to
///     ensure that T is fundamental. This class should be used in lieu
///     of G4THitsCollection<T> when memory is a concern. Atomics are
///     thread-safe and *generally* faster that mutexes (as long as the
///     STL implementation is lock-free) but the synchronization does
///     not come without a cost. If performance is the primary concern,
///     use G4THitsCollection<T> in thread-local instances.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......




#ifndef G4TAtomicHitsCollection_h
#define G4TAtomicHitsCollection_h 1

#include "G4VHitsCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4atomic.hh"

#include <deque>
#include <type_traits>

// class description:
//
//  This is a template class of hits collection and parametrized by
// The concrete class of G4VHit. This is a uniform collection for
// a particular concrete hit class objects.
//  An intermediate layer class G4HitsCollection appeared in this
// header file is used just for G4Allocator, because G4Allocator
// cannot be instansiated with a template class. Thus G4HitsCollection
// class MUST NOT be directly used by the user.

/*class G4HitsCollection : public G4VHitsCollection
{
  public:
      G4HitsCollection();
      G4HitsCollection(G4String detName,G4String colNam);
      virtual ~G4HitsCollection();
      G4bool operator==(const G4HitsCollection &right) const;

  protected:
      void* theCollection;
};*/


template <class T>
class G4TAtomicHitsCollection : public G4VHitsCollection
{
protected:
  static_assert(std::is_fundamental<T>::value,
                "G4TAtomicHitsCollection must use fundamental type");

public:
  typedef T                                 base_type;
  typedef G4atomic<T>                       value_type;
  typedef typename std::deque<value_type*>  container_type;

public:
  G4TAtomicHitsCollection();

public:
  // with description
  G4TAtomicHitsCollection(G4String detName, G4String colNam);

  // constructor.
public:
  virtual ~G4TAtomicHitsCollection();
  G4bool operator==(const G4TAtomicHitsCollection<T> &right) const;

  //inline void *operator new(size_t);
  //inline void operator delete(void* anHC);

public: // with description
  virtual void DrawAllHits();
  virtual void PrintAllHits();
  //  These two methods invokes Draw() and Print() methods of all of
  // hit objects stored in this collection, respectively.

public: // with description
  inline value_type* operator[](size_t i) const
  {
    return (*theCollection)[i];
  }
  //  Returns a pointer to a concrete hit object.
  inline container_type* GetVector() const
  {
    return theCollection;
  }
  //  Returns a collection vector.
  inline G4int insert(T* aHit)
  {
    G4AutoLock l(&fMutex);
    theCollection->push_back(aHit);
    return theCollection->size();
  }
  //  Insert a hit object. Total number of hit objects stored in this
  // collection is returned.
  inline G4int entries() const
  {
    G4AutoLock l(&fMutex);
    return theCollection->size();
  }
  //  Returns the number of hit objects stored in this collection

public:
  virtual G4VHit* GetHit(size_t i) const
  {
    return (*theCollection)[i];
  }
  virtual size_t GetSize() const
  {
    G4AutoLock l(&fMutex);
    return theCollection->size();
  }

protected:
  container_type* theCollection;
  G4Mutex         fMutex;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <class T>
G4TAtomicHitsCollection<T>::G4TAtomicHitsCollection()
  : theCollection(new container_type)
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <class T>
G4TAtomicHitsCollection<T>::G4TAtomicHitsCollection(G4String detName,
                                                    G4String colNam)
  : G4VHitsCollection(detName,colNam),
    theCollection(new container_type)
{ }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <class T> G4TAtomicHitsCollection<T>::~G4TAtomicHitsCollection()
{
  for(size_t i = 0; i < theCollection->size(); i++)
    delete (*theCollection)[i];
  theCollection->clear();
  delete theCollection;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <class T>
G4bool G4TAtomicHitsCollection<T>
::operator==(const G4TAtomicHitsCollection<T> &right) const
{
    return (collectionName == right.collectionName);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <class T>
void G4TAtomicHitsCollection<T>::DrawAllHits()
{
  G4AutoLock l(&fMutex);
  for(size_t i = 0; i < theCollection->size(); i++)
    (*theCollection)[i]->Draw();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
template <class T>
void G4TAtomicHitsCollection<T>::PrintAllHits()
{
  G4AutoLock l(&fMutex);
  for(size_t i = 0; i < theCollection->size(); i++)
    (*theCollection)[i]->Print();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
