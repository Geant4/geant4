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

#ifndef G4HCofThisEvent_h
#define G4HCofThisEvent_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4VHitsCollection.hh"
#include <vector>

// class description:
//
//  This is a class which stores hits collections generated at one event.
// This class is exclusively constructed by G4SDManager when the first
// hits collection of an event is passed to the manager, and this class
// object is deleted by G4RunManager when a G4Event class object is deleted.
//  Almost all public methods must be used by Geant4 kernel classes and
// the user should not invoke them. The user can use two const methods,
// GetHC() and GetNumberOfCollections() for accessing to the stored hits
// collection(s).

class G4HCofThisEvent
{
 public:
  G4HCofThisEvent();
  G4HCofThisEvent(G4int cap);
  ~G4HCofThisEvent();
  inline void* operator new(size_t);
  inline void operator delete(void* anHCoTE);

  void AddHitsCollection(G4int HCID, G4VHitsCollection* aHC);

  G4HCofThisEvent(const G4HCofThisEvent&);
  G4HCofThisEvent& operator=(const G4HCofThisEvent&);

 private:
  std::vector<G4VHitsCollection*>* HC;

 public:  // with description
  inline G4VHitsCollection* GetHC(G4int i) { return (*HC)[i]; }
  //  Returns a pointer to a hits collection. Null will be returned
  // if the particular collection is not stored at the current event.
  // The integer argument is ID number which is assigned by G4SDManager
  // and the number can be obtained by G4SDManager::GetHitsCollectionID()
  // method.
  inline G4int GetNumberOfCollections()
  {
    G4int n = 0;
    for(size_t i = 0; i < HC->size(); i++)
    {
      if((*HC)[i])
        n++;
    }
    return n;
  }
  //  Returns the number of hits collections which are stored in this class
  // object.
 public:
  inline size_t GetCapacity() { return HC->size(); }
};

#if defined G4DIGI_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4HCofThisEvent>*& anHCoTHAllocator_G4MT_TLS_();
#else
extern G4DLLIMPORT G4Allocator<G4HCofThisEvent>*& anHCoTHAllocator_G4MT_TLS_();
#endif

inline void* G4HCofThisEvent::operator new(size_t)
{
  if(anHCoTHAllocator_G4MT_TLS_() == nullptr)
  {
    anHCoTHAllocator_G4MT_TLS_() = new G4Allocator<G4HCofThisEvent>;
  }
  return (void*) anHCoTHAllocator_G4MT_TLS_()->MallocSingle();
}

inline void G4HCofThisEvent::operator delete(void* anHCoTH)
{
  anHCoTHAllocator_G4MT_TLS_()->FreeSingle((G4HCofThisEvent*) anHCoTH);
}

#endif
