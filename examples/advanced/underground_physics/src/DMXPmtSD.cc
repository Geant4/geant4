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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// PmtSD (sensitive PMT) program
// --------------------------------------------------------------

#include "DMXPmtSD.hh"

#include "DMXDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"


DMXPmtSD::DMXPmtSD(G4String name, DMXDetectorConstruction* DMXSD) 
  :G4VSensitiveDetector(name),DMXDetector(DMXSD) {

  G4String HCname="pmtCollection";
  collectionName.insert(HCname);
}


DMXPmtSD::~DMXPmtSD() {;}


////////////////////////////////////////////////////////////////////////////
void DMXPmtSD::Initialize(G4HCofThisEvent* HCE) {

  pmtCollection = new DMXPmtHitsCollection
    (SensitiveDetectorName,collectionName[0]); 

  HitID = -1;


}



////////////////////////////////////////////////////////////////////////////
G4bool DMXPmtSD::ProcessHits
  (G4Step* aStep, G4TouchableHistory* ROhist){

  // make known hit position
  DMXPmtHit* aPmtHit = new DMXPmtHit();
  aPmtHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  aPmtHit->SetTime(aStep->GetPostStepPoint()->GetGlobalTime());
  HitID = pmtCollection->insert(aPmtHit);

  return true;
 
}



////////////////////////////////////////////////////////////////////////////
void DMXPmtSD::EndOfEvent(G4HCofThisEvent* HCE) {

  G4String HCname = collectionName[0];

  static G4int HCID = -1;
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID,pmtCollection);
  
  G4int nHits = pmtCollection->entries();
  if (verboseLevel>=1) {
    G4cout << "     PMT collection: " << nHits << " hits" << G4endl;
    if (verboseLevel>=2)
      pmtCollection->PrintAllHits();
  }


}


////////////////////////////////////////////////////////////////////////////
void DMXPmtSD::clear()    {;}


void DMXPmtSD::DrawAll()  {;}


void DMXPmtSD::PrintAll() {;}




