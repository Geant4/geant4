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
// $Id: ParMarshaledObj.cc,v 1.1 2002-03-05 15:22:03 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------
// When we implement track level parallelism, we will use G4PrimaryTransformer.
// We will also need to merge hits from several slaves. To do so, we will need
// the following:
// - MarshaledHCofThisEvent must also have a private const member:
//     const G4bool canMergeHits = false;
//   and a virtual method (not pure)
//     virtual void insertOrMerge(ExN02TrackerHit *);
// - Then can create a derived class with constructor:
//     inline DerivedClassConstructor()
//       : canMergeHits(true),MarshaledHCofThisEventWithTemplate() {}
// - And insertOrMerge() should be defined in the derived class
//   The logic for insertOrMerge() can usually be copied from
//   (for example) ExN02TrackerSD::ProcessHits
// --------------------------------------------------------------------

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4ios.hh"

#include "ParMarshaledObj.hh"

MarshaledHCofThisEvent::MarshaledHCofThisEvent(G4Event *anEvent)
{ MarshalHitsCollection(anEvent); }

MarshaledHCofThisEvent::MarshaledHCofThisEvent()
{ MarshalHitsCollection( G4RunManager::GetRunManager()-> GetCurrentEvent() ); }

void MarshaledHCofThisEvent::MarshalHitsCollection(const G4Event *anEvent)
{
// MarshaledHCofThisEvent() can be in class ParMarshaledObj
// It can then be dispatched based on HCname to submarshaling
//   for correct class.

  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4int numHitsColl = HCE->GetNumberOfCollections();
  Marshal( numHitsColl );
  for (int i = 0; i < numHitsColl; i++) {
    // Use base class to find HCname.
    G4VHitsCollection *aBaseHC = HCE->GetHC(i);
    G4String SDname = aBaseHC->GetSDname();
    G4String HCname = aBaseHC->GetName();
    Marshal( SDname );
    Marshal( HCname );
#ifdef TOPC_DEBUG
    G4cout << "Marshaled HCname: " << HCname << G4endl;
    G4cout << "Marshaled SDname: " << SDname << G4endl;
#endif

    MarshalHitsCollection(aBaseHC, HCname);
  }
}

// NOTES TO MYSELF:
// G4SDManager::GetSDMpointer()->FindSensitiveDetector(SDname);

// IMPORTANT:  all calls to Unmarshal must be in same order as
//             the calls to Marshal in the constructor.
void MarshaledHCofThisEvent::UnmarshalSlaveHCofThisEvent()
{
   /*
   G4HCofThisEvent *HCE = G4RunManager::GetRunManager()->
			  GetCurrentEvent()->
		          GetHCofThisEvent();
   */
   G4int numHitsColl;
   Unmarshal( numHitsColl );
   for (int i = 0; i < numHitsColl; i++) {
     G4String SDname;
     G4String HCname;
     Unmarshal( SDname );
     Unmarshal( HCname );
#ifdef TOPC_DEBUG
     G4cout << "HCname: " << HCname << G4endl;
     G4cout << "SDname: " << SDname << G4endl;
#endif

     UnmarshalSlaveHitsCollection(HCname, SDname);
  }
}
