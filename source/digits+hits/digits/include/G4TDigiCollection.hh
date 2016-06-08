//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TDigiCollection.hh,v 1.8.2.1 2001/06/28 19:07:50 gunter Exp $
// GEANT4 tag $Name:  $
//

#ifndef G4TDigiCollection_h
#define G4TDigiCollection_h 1

#include "G4VDigiCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
//#include "g4rw/tpordvec.h"
#include "g4std/vector"

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
      { return (*((G4std::vector<T*>*)theCollection))[i]; }
      //  Returns a pointer to a concrete digi object.
      inline G4std::vector<T*>* GetVector() const
      { return (G4std::vector<T*>*)theCollection; }
      //  Returns a collection vector.
      inline int insert(T* aHit)
      {
        G4std::vector<T*>*theDigiCollection 
          = (G4std::vector<T*>*)theCollection;
        theDigiCollection->push_back(aHit);
        return theDigiCollection->size();
      }
      //  Insert a digi object. Total number of digi objects stored in this
      // collection is returned.
      inline int entries() const
      {
        G4std::vector<T*>*theDigiCollection 
          = (G4std::vector<T*>*)theCollection;
        return theDigiCollection->size();
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
  G4std::vector<T*> * theDigiCollection
    = new G4std::vector<T*>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::G4TDigiCollection(G4String detName,G4String colNam)
: G4DigiCollection(detName,colNam)
{ 
  G4std::vector<T*> * theDigiCollection
    = new G4std::vector<T*>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::~G4TDigiCollection()
{
  G4std::vector<T*> * theDigiCollection 
    = (G4std::vector<T*>*)theCollection;
  //theDigiCollection->clearAndDestroy();
  for(int i=0;i<theDigiCollection->size();i++)
  { delete (*theDigiCollection)[i]; }
  theDigiCollection->clear();
  delete theDigiCollection;
}

template <class T> int G4TDigiCollection<T>::operator==(const G4TDigiCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4TDigiCollection<T>::DrawAllDigi() 
{
  G4std::vector<T*> * theDigiCollection 
    = (G4std::vector<T*>*)theCollection;
  int n = theDigiCollection->size();
  for(int i=0;i<n;i++)
  { (*theDigiCollection)[i]->Draw(); }
}

template <class T> void G4TDigiCollection<T>::PrintAllDigi() 
{
  G4std::vector<T*> * theDigiCollection 
    = (G4std::vector<T*>*)theCollection;
  int n = theDigiCollection->size();
  for(int i=0;i<n;i++)
  { (*theDigiCollection)[i]->Print(); }
}

#endif

