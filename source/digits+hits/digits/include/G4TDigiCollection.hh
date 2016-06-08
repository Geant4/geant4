// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TDigiCollection.hh,v 1.5.2.1.2.1 1999/12/07 20:47:46 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#ifndef G4TDigiCollection_h
#define G4TDigiCollection_h 1

#include "G4VDigiCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
#include "g4rw/tpordvec.h"

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
  public: // with description
      G4TDigiCollection(G4String detName,G4String colNam);
      // Constructor.
  public:
      virtual ~G4TDigiCollection();
      int operator==(const G4TDigiCollection &right) const;
      
      inline void *operator new(size_t);
      inline void operator delete(void* aDC);
  public: // with description
      virtual void DrawAllDigi();
      virtual void PrintAllDigi();
      //  These two methods invokes Draw() and Print() methods of all of
      // digit objects stored in this collection, respectively.

  public: // with description
      inline T* operator[](size_t i) const
      { return (*((G4RWTPtrOrderedVector<T>*)theCollection))[i]; }
      //  Returns a pointer to a concrete digi object.
      inline G4RWTPtrOrderedVector<T>* GetVector() const
      { return (G4RWTPtrOrderedVector<T>*)theCollection; }
      //  Returns a collection vector.
      inline int insert(T* aHit)
      {
        G4RWTPtrOrderedVector<T>*theDigiCollection 
          = (G4RWTPtrOrderedVector<T>*)theCollection;
        theDigiCollection->insert(aHit);
        return theDigiCollection->entries();
      }
      //  Insert a digi object. Total number of digi objects stored in this
      // collection is returned.
      inline int entries() const
      {
        G4RWTPtrOrderedVector<T>*theDigiCollection 
          = (G4RWTPtrOrderedVector<T>*)theCollection;
        return theDigiCollection->entries();
      }
      //  Returns the number of digi objcets stored in this collection.
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
  G4RWTPtrOrderedVector<T> * theDigiCollection
    = new G4RWTPtrOrderedVector<T>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::G4TDigiCollection(G4String detName,G4String colNam)
: G4DigiCollection(detName,colNam)
{ 
  G4RWTPtrOrderedVector<T> * theDigiCollection
    = new G4RWTPtrOrderedVector<T>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::~G4TDigiCollection()
{
  G4RWTPtrOrderedVector<T> * theDigiCollection 
    = (G4RWTPtrOrderedVector<T>*)theCollection;
  theDigiCollection->clearAndDestroy();
  delete theDigiCollection;
}

template <class T> int G4TDigiCollection<T>::operator==(const G4TDigiCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4TDigiCollection<T>::DrawAllDigi() 
{
  G4RWTPtrOrderedVector<T> * theDigiCollection 
    = (G4RWTPtrOrderedVector<T>*)theCollection;
  int n = theDigiCollection->entries();
  for(int i=0;i<n;i++)
  { (*theDigiCollection)[i]->Draw(); }
}

template <class T> void G4TDigiCollection<T>::PrintAllDigi() 
{
  G4RWTPtrOrderedVector<T> * theDigiCollection 
    = (G4RWTPtrOrderedVector<T>*)theCollection;
  int n = theDigiCollection->entries();
  for(int i=0;i<n;i++)
  { (*theDigiCollection)[i]->Print(); }
}

#endif

