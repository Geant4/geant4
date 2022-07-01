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

#ifndef G4DCofThisEvent_h
#define G4DCofThisEvent_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4VDigiCollection.hh"
//#include "g4rw/tpordvec.h"
#include <vector>

// class description:
//
//  This is a class which stores digi collections generated at one event.
// This class is exclusively constructed by G4DigiManager when the first
// digi collection of an event is passed to the manager, and this class
// object is deleted by G4RunManager when a G4Event class object is deleted.
//  Almost all public methods must be used by Geant4 kernel classes and
// the user should not invoke them. The user can use two const methods,
// GetDC() and GetNumberOfCollections() for accessing to the stored digi
// collection(s).

class G4DCofThisEvent
{
 public:
  G4DCofThisEvent();
  G4DCofThisEvent(G4int cap);
  ~G4DCofThisEvent();
  inline void* operator new(size_t);
  inline void operator delete(void* anDCoTE);

  void AddDigiCollection(G4int DCID, G4VDigiCollection* aDC);

  G4DCofThisEvent(const G4DCofThisEvent&);
  G4DCofThisEvent& operator=(const G4DCofThisEvent&);

 private:
  std::vector<G4VDigiCollection*>* DC;

 public:  // with description
  inline G4VDigiCollection* GetDC(G4int i) const { return (*DC)[i]; }
  //  Returns a pointer to a digi collection. Null will be returned
  // if the particular collection is not stored at the current event.
  // The integer argument is ID number which is assigned by G4DigiManager
  // and the number can be obtained by G4DigiManager::GetDigiCollectionID()
  // method.
  inline G4int GetNumberOfCollections() const
  {
    G4int n = 0;
    for(auto dc : *DC)
      if(dc)
        n++;
    return n;
  }
  //  Returns the number of digi collections which are stored in this class
  // object.
 public:
  inline size_t GetCapacity() const { return DC->size(); }
};

#if defined G4DIGI_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4DCofThisEvent>*& anDCoTHAllocator_G4MT_TLS_();
#else
extern G4DLLIMPORT G4Allocator<G4DCofThisEvent>*& anDCoTHAllocator_G4MT_TLS_();
#endif

inline void* G4DCofThisEvent::operator new(size_t)
{
  if(anDCoTHAllocator_G4MT_TLS_() == nullptr)
  {
    anDCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4DCofThisEvent>;
  }
  return (void*) anDCoTHAllocator_G4MT_TLS_()->MallocSingle();
}

inline void G4DCofThisEvent::operator delete(void* anDCoTH)
{
  anDCoTHAllocator_G4MT_TLS_()->FreeSingle((G4DCofThisEvent*) anDCoTH);
}

#endif
