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
/// \file G4THitsCollection.hh
/// \brief Definition of the G4THitsCollection class
//
//
//

#ifndef G4THitsCollection_h
#define G4THitsCollection_h 1

#include "G4VHitsCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"
//#include "g4rw/tpordvec.h"
#include <vector>

//MSH_include_begin
#include "ExN04CalorimeterHit.hh"
#include "ExN04MuonHit.hh"
#include "ExN04TrackerHit.hh"
#include "MarshaledExN04CalorimeterHit.h"
#include "MarshaledExN04MuonHit.h"
#include "MarshaledExN04TrackerHit.h"
//MSH_include_end


// class description:
//
//  This is a template class of hits collection and parametrized by
// The concrete class of G4VHit. This is a uniform collection for
// a particular concrete hit class objects.
//  An intermediate layer class G4HitsCollection appeared in this
// header file is used just for G4Allocator, because G4Allocator
// cannot be instansiated with a template class. Thus G4HitsCollection
// class MUST NOT be directly used by the user.

//MSH_BEGIN
class G4HitsCollection : public G4VHitsCollection
{
  public:
      G4HitsCollection();
      G4HitsCollection(G4String detName,G4String colNam);
      virtual ~G4HitsCollection();
      G4bool operator==(const G4HitsCollection &right) const;

  protected:
      void* theCollection; /*MSH: ptr_as_array
    [elementType: (dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>($THIS)!=NULL) => ExN04CalorimeterHit* 
        | (dynamic_cast<G4THitsCollection<ExN04MuonHit>*>($THIS)!=NULL) => ExN04MuonHit* 
        | true => ExN04TrackerHit*]
    [elementCount: { 
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>($THIS)!=NULL)
        $ELE_COUNT = ((G4THitsCollection<ExN04CalorimeterHit>*)$THIS)->entries(); 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>($THIS)!=NULL)
        $ELE_COUNT = ((G4THitsCollection<ExN04MuonHit>*)$THIS)->entries(); 
     else
        $ELE_COUNT = ((G4THitsCollection<ExN04TrackerHit>*)$THIS)->entries(); 
    }]
    [elementGet: { 
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>($THIS)!=NULL)
        $ELEMENT = (*((G4THitsCollection<ExN04CalorimeterHit>*)$THIS))[$ELE_INDEX]; 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>($THIS)!=NULL)
        $ELEMENT = (*((G4THitsCollection<ExN04MuonHit>*)$THIS))[$ELE_INDEX]; 
     else
        $ELEMENT = (*((G4THitsCollection<ExN04TrackerHit>*)$THIS))[$ELE_INDEX]; 
    }]
    [elementSet: {
     if(dynamic_cast<G4THitsCollection<ExN04CalorimeterHit>*>($THIS)!=NULL)
        ((G4THitsCollection<ExN04CalorimeterHit>*)$THIS)->insert((ExN04CalorimeterHit*)$ELEMENT); 
     else if(dynamic_cast<G4THitsCollection<ExN04MuonHit>*>($THIS)!=NULL)
        ((G4THitsCollection<ExN04MuonHit>*)$THIS)->insert((ExN04MuonHit*)$ELEMENT); 
     else
        ((G4THitsCollection<ExN04TrackerHit>*)$THIS)->insert((ExN04TrackerHit*)$ELEMENT); 
    }]	*/

};
//MSH_END


#if defined G4DIGI_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4HitsCollection> anHCAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4HitsCollection> anHCAllocator;
#endif

//MSH_BEGIN
template <class T> class G4THitsCollection : public G4HitsCollection 
{
  public:
      G4THitsCollection();
  public: // with description
      G4THitsCollection(G4String detName,G4String colNam);
      // constructor.
  public:
      virtual ~G4THitsCollection();
      G4bool operator==(const G4THitsCollection<T> &right) const;
      
      inline void *operator new(size_t);
      inline void operator delete(void* anHC);
  public: // with description
      virtual void DrawAllHits();
      virtual void PrintAllHits();
      //  These two methods invokes Draw() and Print() methods of all of
      // hit objects stored in this collection, respectively.

  public: // with description
      inline T* operator[](size_t i) const
      { return (*((std::vector<T*>*)theCollection))[i]; }
      //  Returns a pointer to a concrete hit object.
      inline std::vector<T*>* GetVector() const
      { return (std::vector<T*>*)theCollection; }
      //  Returns a collection vector.
      inline G4int insert(T* aHit)
      {
        std::vector<T*>*theHitsCollection 
          = (std::vector<T*>*)theCollection;
        theHitsCollection->push_back(aHit);
        return theHitsCollection->size();
      }
      //  Insert a hit object. Total number of hit objects stored in this
      // collection is returned.
      inline G4int entries() const
      {
        std::vector<T*>*theHitsCollection
          = (std::vector<T*>*)theCollection;
        return theHitsCollection->size();
      }
      //  Returns the number of hit objects stored in this collection

  public:
      virtual G4VHit* GetHit(size_t i) const
      { return (*((std::vector<T*>*)theCollection))[i]; }
      virtual size_t GetSize() const
      { return ((std::vector<T*>*)theCollection)->size(); }

  // MSH_superclass : G4HitsCollection

};
//MSH_END

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
  std::vector<T*> * theHitsCollection
    = new std::vector<T*>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::G4THitsCollection(G4String detName,G4String colNam)
: G4HitsCollection(detName,colNam)
{ 
  std::vector<T*> * theHitsCollection
    = new std::vector<T*>;
  theCollection = (void*)theHitsCollection;
}

template <class T> G4THitsCollection<T>::~G4THitsCollection()
{
  std::vector<T*> * theHitsCollection 
    = (std::vector<T*>*)theCollection;
  //theHitsCollection->clearAndDestroy();
  for(size_t i=0;i<theHitsCollection->size();i++)
  { delete (*theHitsCollection)[i]; }
  theHitsCollection->clear();
  delete theHitsCollection;
}

template <class T> G4bool G4THitsCollection<T>::operator==(const G4THitsCollection<T> &right) const
{ return (collectionName==right.collectionName); }

template <class T> void G4THitsCollection<T>::DrawAllHits() 
{
  std::vector<T*> * theHitsCollection 
    = (std::vector<T*>*)theCollection;
  size_t n = theHitsCollection->size();
  for(size_t i=0;i<n;i++)
  { (*theHitsCollection)[i]->Draw(); }
}

template <class T> void G4THitsCollection<T>::PrintAllHits() 
{
  std::vector<T*> * theHitsCollection 
    = (std::vector<T*>*)theCollection;
  size_t n = theHitsCollection->size();
  for(size_t i=0;i<n;i++)
  { (*theHitsCollection)[i]->Print(); }
}

#endif

