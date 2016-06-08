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
// $Id: G4THitsCollection.hh,v 1.7.2.1 2001/06/28 19:07:51 gunter Exp $
// GEANT4 tag $Name:  $
//

#ifndef G4THitsCollection_h
#define G4THitsCollection_h 1

#include "G4VHitsCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
//#include "g4rw/tpordvec.h"
#include "g4std/vector"

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
      { return (*((G4std::vector<T*>*)theCollection))[i]; }
      //  Returns a pointer to a concrete hit object.
      inline G4std::vector<T*>* GetVector() const
      { return (G4std::vector<T*>*)theCollection; }
      //  Returns a collection vector.
      inline int insert(T* aHit)
      {
        G4std::vector<T*>*theHitsCollection 
          = (G4std::vector<T*>*)theCollection;
        theHitsCollection->push_back(aHit);
        return theHitsCollection->size();
      }
      //  Insert a hit object. Total number of hit objects stored in this
      // collection is returned.
      inline int entries() const
      {
        G4std::vector<T*>*theHitsCollection
          = (G4std::vector<T*>*)theCollection;
        return theHitsCollection->size();
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
  G4std::vector<T*> * theHitsCollection
    = new G4std::vector<T*>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::G4THitsCollection(G4String detName,G4String colNam)
: G4HitsCollection(detName,colNam)
{ 
  G4std::vector<T*> * theHitsCollection
    = new G4std::vector<T*>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::~G4THitsCollection()
{
  G4std::vector<T*> * theHitsCollection 
    = (G4std::vector<T*>*)theCollection;
  //theHitsCollection->clearAndDestroy();
  for(int i=0;i<theHitsCollection->size();i++)
  { delete (*theHitsCollection)[i]; }
  theHitsCollection->clear();
  delete theHitsCollection;
}

template <class T> int G4THitsCollection<T>::operator==(const G4THitsCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4THitsCollection<T>::DrawAllHits() 
{
  G4std::vector<T*> * theHitsCollection 
    = (G4std::vector<T*>*)theCollection;
  int n = theHitsCollection->size();
  for(int i=0;i<n;i++)
  { (*theHitsCollection)[i]->Draw(); }
}

template <class T> void G4THitsCollection<T>::PrintAllHits() 
{
  G4std::vector<T*> * theHitsCollection 
    = (G4std::vector<T*>*)theCollection;
  int n = theHitsCollection->size();
  for(int i=0;i<n;i++)
  { (*theHitsCollection)[i]->Print(); }
}

#endif

