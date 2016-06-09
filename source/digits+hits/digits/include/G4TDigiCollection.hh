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
// $Id: G4TDigiCollection.hh,v 1.11 2003/06/16 16:50:20 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
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
      G4DigiCollection(G4String detName,G4String colNam);
      virtual ~G4DigiCollection();
      G4int operator==(const G4DigiCollection &right) const;

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
      G4int operator==(const G4TDigiCollection &right) const;
      
      inline void *operator new(size_t);
      inline void operator delete(void* aDC);
  public: // with description
      virtual void DrawAllDigi();
      virtual void PrintAllDigi();
      //  These two methods invokes Draw() and Print() methods of all of
      // digit objects stored in this collection, respectively.

  public: // with description
      inline T* operator[](size_t i) const
      { return (*((std::vector<T*>*)theCollection))[i]; }
      //  Returns a pointer to a concrete digi object.
      inline std::vector<T*>* GetVector() const
      { return (std::vector<T*>*)theCollection; }
      //  Returns a collection vector.
      inline G4int insert(T* aHit)
      {
        std::vector<T*>*theDigiCollection 
          = (std::vector<T*>*)theCollection;
        theDigiCollection->push_back(aHit);
        return theDigiCollection->size();
      }
      //  Insert a digi object. Total number of digi objects stored in this
      // collection is returned.
      inline G4int entries() const
      {
        std::vector<T*>*theDigiCollection 
          = (std::vector<T*>*)theCollection;
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
  std::vector<T*> * theDigiCollection
    = new std::vector<T*>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::G4TDigiCollection(G4String detName,G4String colNam)
: G4DigiCollection(detName,colNam)
{ 
  std::vector<T*> * theDigiCollection
    = new std::vector<T*>;
  theCollection = (void*)theDigiCollection;
}

template <class T> G4TDigiCollection<T>::~G4TDigiCollection()
{
  std::vector<T*> * theDigiCollection 
    = (std::vector<T*>*)theCollection;
  //theDigiCollection->clearAndDestroy();
  for(size_t i=0;i<theDigiCollection->size();i++)
  { delete (*theDigiCollection)[i]; }
  theDigiCollection->clear();
  delete theDigiCollection;
}

template <class T> G4int G4TDigiCollection<T>::operator==(const G4TDigiCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4TDigiCollection<T>::DrawAllDigi() 
{
  std::vector<T*> * theDigiCollection 
    = (std::vector<T*>*)theCollection;
  size_t n = theDigiCollection->size();
  for(size_t i=0;i<n;i++)
  { (*theDigiCollection)[i]->Draw(); }
}

template <class T> void G4TDigiCollection<T>::PrintAllDigi() 
{
  std::vector<T*> * theDigiCollection 
    = (std::vector<T*>*)theCollection;
  size_t n = theDigiCollection->size();
  for(size_t i=0;i<n;i++)
  { (*theDigiCollection)[i]->Print(); }
}

#endif

