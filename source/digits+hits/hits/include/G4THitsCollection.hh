// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4THitsCollection.hh,v 2.1 1998/07/12 02:53:42 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4THitsCollection_h
#define G4THitsCollection_h 1

#include "G4VHitsCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
#include <rw/tpordvec.h>

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
      G4THitsCollection(G4String detName,G4String colNam);
      virtual ~G4THitsCollection();
      int operator==(const G4THitsCollection<T> &right) const;
      
      inline void *operator new(size_t);
      inline void operator delete(void* anHC);

      virtual void DrawAllHits();
      virtual void PrintAllHits();

  public:
      inline T* operator[](size_t i)
      { return (*((RWTPtrOrderedVector<T>*)theCollection))[i]; }
      inline RWTPtrOrderedVector<T>* GetVector()
      { return (RWTPtrOrderedVector<T>*)theCollection; }
      inline int insert(T* aHit)
      {
        RWTPtrOrderedVector<T>*theHitsCollection 
          = (RWTPtrOrderedVector<T>*)theCollection;
        theHitsCollection->insert(aHit);
        return theHitsCollection->entries();
      }
      inline int entries()
      {
        RWTPtrOrderedVector<T>*theHitsCollection
          = (RWTPtrOrderedVector<T>*)theCollection;
        return theHitsCollection->entries();
      }

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
  RWTPtrOrderedVector<T> * theHitsCollection
    = new RWTPtrOrderedVector<T>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::G4THitsCollection(G4String detName,G4String colNam)
: G4HitsCollection(detName,colNam)
{ 
  RWTPtrOrderedVector<T> * theHitsCollection
    = new RWTPtrOrderedVector<T>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::~G4THitsCollection()
{
  RWTPtrOrderedVector<T> * theHitsCollection 
    = (RWTPtrOrderedVector<T>*)theCollection;
  theHitsCollection->clearAndDestroy();
  delete theHitsCollection;
}

template <class T> int G4THitsCollection<T>::operator==(const G4THitsCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4THitsCollection<T>::DrawAllHits()
{
  RWTPtrOrderedVector<T> * theHitsCollection 
    = (RWTPtrOrderedVector<T>*)theCollection;
  int n = theHitsCollection->entries();
  for(int i=0;i<n;i++)
  { (*theHitsCollection)[i]->Draw(); }
}

template <class T> void G4THitsCollection<T>::PrintAllHits()
{
  RWTPtrOrderedVector<T> * theHitsCollection 
    = (RWTPtrOrderedVector<T>*)theCollection;
  int n = theHitsCollection->entries();
  for(int i=0;i<n;i++)
  { (*theHitsCollection)[i]->Print(); }
}

#endif

