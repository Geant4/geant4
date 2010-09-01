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

#include "Test2PhantomSD.hh"
#include "Test2PhantomHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"

Test2PhantomSD::Test2PhantomSD(G4String name)
  :G4VSensitiveDetector(name) {

  G4String HCname;
  collectionName.insert(HCname = "PhantomCollection");

}

Test2PhantomSD::~Test2PhantomSD() {
  ;
}

void Test2PhantomSD::Initialize(G4HCofThisEvent *) {

  fPhantomCollection = new Test2PhantomHitsCollection(SensitiveDetectorName,
						      collectionName[0]); 
  verboseLevel = 0;
}

G4bool Test2PhantomSD::ProcessHits(G4Step * aStep, G4TouchableHistory *) {

  // total energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double trklen = aStep->GetStepLength();
  if(verboseLevel > 1) G4cout << "Next step edep [MeV] = " << edep/MeV << G4endl;
  if(edep > 0. || trklen > 0.) {

    G4TouchableHistory * hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4int copyIDinX = hist->GetReplicaNumber(0);
    G4int copyIDinY = hist->GetReplicaNumber(1);
    G4int copyIDinZ = hist->GetReplicaNumber(2);

    Test2PhantomHit* phantomHit
	= new Test2PhantomHit(copyIDinX, copyIDinY, copyIDinZ);
    phantomHit->SetEdep(edep);
    phantomHit->SetTrackLength(trklen);
    phantomHit->SetParticleName(aStep->GetTrack()->GetParticleDefinition()->GetParticleName());
    fPhantomCollection->insert(phantomHit);
    if(verboseLevel > 0) {
      G4cout << " A hit of Test2PhantomHit is created in copy-id (" 
	     << copyIDinX << ", " << copyIDinY << ", " << copyIDinZ
	     << ")" << G4endl;
    }
  }

  return true;
}

void Test2PhantomSD::EndOfEvent(G4HCofThisEvent * HCE) {

  static G4int HCID = -1;
  if(HCID < 0) { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, fPhantomCollection );
}

void Test2PhantomSD::clear() {
  ;
} 

void Test2PhantomSD::DrawAll() {
  ;
} 

void Test2PhantomSD::PrintAll() {
  ;
} 

