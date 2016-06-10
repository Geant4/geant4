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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
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


DMXPmtSD::DMXPmtSD(G4String name) 
  :G4VSensitiveDetector(name) {

  G4String HCname="pmtCollection";
  collectionName.insert(HCname);
}


DMXPmtSD::~DMXPmtSD() {;}


////////////////////////////////////////////////////////////////////////////
void DMXPmtSD::Initialize(G4HCofThisEvent*) {

  pmtCollection = new DMXPmtHitsCollection
    (SensitiveDetectorName,collectionName[0]); 

  HitID = -1;


}



////////////////////////////////////////////////////////////////////////////
G4bool DMXPmtSD::ProcessHits
  (G4Step* aStep, G4TouchableHistory*){

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




