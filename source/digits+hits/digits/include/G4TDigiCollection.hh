// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TDigiCollection.hh,v 1.1 1999/01/07 16:06:28 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

#ifndef G4TDigiCollection_h
#define G4TDigiCollection_h 1

#include "G4VDigiCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
#include <rw/tpordvec.h>

class G4DigiCollection : public G4VDigiCollection
{
  public:
      G4DigiCollection();
      G4DigiCollection(G4String detName,G4String colNam);
      virtual ~G4DigiCollection();
      int operator==(const G4DigiCollection &right) const;

  protected:
      void* theCollection;
};

extern G4Allocator<G4DigiCollection> aDCAllocator;

template <class T> class G4TDigiCollection : public G4DigiCollection 
{
  public:
      G4TDigiCollection();
      G4TDigiCollection(G4String detName,G4String colNam);
      virtual ~G4TDigiCollection();
      int operator==(const G4TDigiCollection &right) const;
      
      inline void *operator new(size_t);
      inline void operator delete(void* aDC);

      virtual void DrawAllDigi();
      virtual void PrintAllDigi();

  public:
      inline T* operator[](size_t i)
      { return (*((RWTPtrOrderedVector<T>*)theCollection))[i]; }
      inline RWTPtrOrderedVector<T>* GetVector()
      { return (RWTPtrOrderedVector<T>*)theCollection; }
      inline int insert(T* aHit)
      {
        RWTPtrOrderedVector<T>*theDigiCollection 
          = (RWTPtrOrderedVector<T>*)theCollection;
        theDigiCollection->insert(aHit);
        return theDigiCollection->entries();
      }
};

template <class T> inline void* G4TDigiCollection<T>::operator new(size_t)
{
  void* aDC;
  aDC = (void*)aDCAllocator.MallocSingle();
  return aDC;
}

template <class T> inline void G4TDigiCollection<T>::operator delete(void* aDC)
{
  aDCAllocator.FreeSingle((G4DigiCollection*)aDC);
}

template <class T> G4TDigiCollection<T>::G4TDigiCollection()
{ 
  RWTPtrOrderedVector<T> * theDigiCollection
    = new RWTPtrOrderedVector<T>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::G4TDigiCollection(G4String detName,G4String colNam)
: G4DigiCollection(detName,colNam)
{ 
  RWTPtrOrderedVector<T> * theDigiCollection
    = new RWTPtrOrderedVector<T>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::~G4TDigiCollection()
{
  RWTPtrOrderedVector<T> * theDigiCollection 
    = (RWTPtrOrderedVector<T>*)theCollection;
  theDigiCollection->clearAndDestroy();
  delete theDigiCollection;
}

template <class T> int G4TDigiCollection<T>::operator==(const G4TDigiCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4TDigiCollection<T>::DrawAllDigi()
{
  RWTPtrOrderedVector<T> * theDigiCollection 
    = (RWTPtrOrderedVector<T>*)theCollection;
  int n = theDigiCollection->entries();
  for(int i=0;i<n;i++)
  { (*theDigiCollection)[i]->Draw(); }
}

template <class T> void G4TDigiCollection<T>::PrintAllDigi()
{
  RWTPtrOrderedVector<T> * theDigiCollection 
    = (RWTPtrOrderedVector<T>*)theCollection;
  int n = theDigiCollection->entries();
  for(int i=0;i<n;i++)
  { (*theDigiCollection)[i]->Print(); }
}

#endif

