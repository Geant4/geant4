// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4THitsCollection.hh,v 1.4.2.1.2.1 1999/12/07 20:47:47 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#ifndef G4THitsCollection_h
#define G4THitsCollection_h 1

#include "G4VHitsCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
#include "g4rw/tpordvec.h"

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
      G4HitsCollection();
      G4HitsCollection(G4String detName,G4String colNam);
      virtual ~G4HitsCollection();
      int operator==(const G4HitsCollection &right) const;

  protected:
      void* theCollection;
};

extern G4Allocator<G4HitsCollection> anHCAllocator;

template <class T> class G4THitsCollection : public G4HitsCollection 
{
  public:
      G4THitsCollection();
  public: // with description
      G4THitsCollection(G4String detName,G4String colNam);
      // constructor.
  public:
      virtual ~G4THitsCollection();
      int operator==(const G4THitsCollection<T> &right) const;
      
      inline void *operator new(size_t);
      inline void operator delete(void* anHC);
  public: // with description
      virtual void DrawAllHits();
      virtual void PrintAllHits();
      //  These two methods invokes Draw() and Print() methods of all of
      // hit objects stored in this collection, respectively.

  public: // with description
      inline T* operator[](size_t i) const
      { return (*((G4RWTPtrOrderedVector<T>*)theCollection))[i]; }
      //  Returns a pointer to a concrete hit object.
      inline G4RWTPtrOrderedVector<T>* GetVector() const
      { return (G4RWTPtrOrderedVector<T>*)theCollection; }
      //  Returns a collection vector.
      inline int insert(T* aHit)
      {
        G4RWTPtrOrderedVector<T>*theHitsCollection 
          = (G4RWTPtrOrderedVector<T>*)theCollection;
        theHitsCollection->insert(aHit);
        return theHitsCollection->entries();
      }
      //  Insert a hit object. Total number of hit objects stored in this
      // collection is returned.
      inline int entries() const
      {
        G4RWTPtrOrderedVector<T>*theHitsCollection
          = (G4RWTPtrOrderedVector<T>*)theCollection;
        return theHitsCollection->entries();
      }
      //  Returns the number of hit objects stored in this collection

};

template <class T> inline void* G4THitsCollection<T>::operator new(size_t)
{
  void* anHC;
  anHC = (void*)anHCAllocator.MallocSingle();
  return anHC;
}

template <class T> inline void G4THitsCollection<T>::operator delete(void* anHC)
{
  anHCAllocator.FreeSingle((G4HitsCollection*)anHC);
}

template <class T> G4THitsCollection<T>::G4THitsCollection()
{ 
  G4RWTPtrOrderedVector<T> * theHitsCollection
    = new G4RWTPtrOrderedVector<T>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::G4THitsCollection(G4String detName,G4String colNam)
: G4HitsCollection(detName,colNam)
{ 
  G4RWTPtrOrderedVector<T> * theHitsCollection
    = new G4RWTPtrOrderedVector<T>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::~G4THitsCollection()
{
  G4RWTPtrOrderedVector<T> * theHitsCollection 
    = (G4RWTPtrOrderedVector<T>*)theCollection;
  theHitsCollection->clearAndDestroy();
  delete theHitsCollection;
}

template <class T> int G4THitsCollection<T>::operator==(const G4THitsCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4THitsCollection<T>::DrawAllHits() 
{
  G4RWTPtrOrderedVector<T> * theHitsCollection 
    = (G4RWTPtrOrderedVector<T>*)theCollection;
  int n = theHitsCollection->entries();
  for(int i=0;i<n;i++)
  { (*theHitsCollection)[i]->Draw(); }
}

template <class T> void G4THitsCollection<T>::PrintAllHits() 
{
  G4RWTPtrOrderedVector<T> * theHitsCollection 
    = (G4RWTPtrOrderedVector<T>*)theCollection;
  int n = theHitsCollection->entries();
  for(int i=0;i<n;i++)
  { (*theHitsCollection)[i]->Print(); }
}

#endif

