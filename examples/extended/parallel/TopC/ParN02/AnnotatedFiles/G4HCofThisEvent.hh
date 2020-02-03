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
/// \file G4HCofThisEvent.hh
/// \brief Definition of the G4HCofThisEvent class
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


//vietha 2003/04/26 , include files needed to copy to the beginning of MarshaledG4HCofThisEvent.h
// in order to compile
//MSH_include_begin
#include "G4SDManager.hh"
#include "G4THitsCollection.hh"
#include "ExN02TrackerHit.hh"
#include "MarshaledExN02TrackerHit.h"
#include "MarshaledG4THitsCollection.h"
#include "MarshaledG4VHitsCollection.h"
//MSH_include_end
//MSH_BEGIN
class G4HCofThisEvent 
{
  public:
      G4HCofThisEvent();
      G4HCofThisEvent(G4int cap);
      ~G4HCofThisEvent();
      inline void *operator new(size_t);
      inline void operator delete(void* anHCoTE);

      void AddHitsCollection(G4int HCID,G4VHitsCollection * aHC);

  private:
      std::vector<G4VHitsCollection*> * HC; /*MSH: ptr_as_array
  [elementType: G4VHitsCollection*]
  [elementCount: { $ELE_COUNT = $THIS->GetNumberOfCollections(); }]
  [elementGet: { $ELEMENT = $THIS->GetHC($ELE_INDEX); }]
  [elementSet: { $THIS->AddHitsCollection($ELE_INDEX, $ELEMENT); }]  */

  public: // with description
      inline G4VHitsCollection* GetHC(G4int i)
      { return (*HC)[i]; }
      //  Returns a pointer to a hits collection. Null will be returned
      // if the particular collection is not stored at the current event.
      // The integer argument is ID number which is assigned by G4SDManager
      // and the number can be obtained by G4SDManager::GetHitsCollectionID()
      // method.
      inline G4int GetNumberOfCollections()
      {
        G4int n = 0;
        for(size_t i=0;i<HC->size();i++)
        {
          if((*HC)[i]) n++;
        }
        return n;
      }
      //  Returns the number of hits collections which are stored in this class
      // object.
  public:
      inline G4int GetCapacity()
      {
        return HC->size();
      }
};
//MSH_END

#if defined G4DIGI_ALLOC_EXPORT
  extern G4DLLEXPORT G4Allocator<G4HCofThisEvent> anHCoTHAllocator;
#else
  extern G4DLLIMPORT G4Allocator<G4HCofThisEvent> anHCoTHAllocator;
#endif

inline void* G4HCofThisEvent::operator new(size_t)
{
  void* anHCoTH;
  anHCoTH = (void*)anHCoTHAllocator.MallocSingle();
  return anHCoTH;
}

inline void G4HCofThisEvent::operator delete(void* anHCoTH)
{
  anHCoTHAllocator.FreeSingle((G4HCofThisEvent*)anHCoTH);
}


#endif

