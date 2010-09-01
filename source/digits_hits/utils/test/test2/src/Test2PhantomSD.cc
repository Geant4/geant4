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
  :G4VSensitiveDetector(name),
   fNumberOfCellsInX(10), fNumberOfCellsInY(10), fNumberOfCellsInZ(10) {

  G4String HCname;
  collectionName.insert(HCname = "PhantomCollection");

}

Test2PhantomSD::~Test2PhantomSD() {
  ;
}

void Test2PhantomSD::Initialize(G4HCofThisEvent *) {

  fPhantomCollection = new Test2PhantomHitsCollection(SensitiveDetectorName,
						      collectionName[0]); 
  for(G4int i = 0; i < fNumberOfCellsInZ; i++) {
    for(G4int j = 0; j < fNumberOfCellsInY; j++) {
      for(G4int k = 0; k < fNumberOfCellsInX; k++) {
	fCellID[i][j][k] = -1;
      }
    }
  }
  verboseLevel = 0;
}

G4bool Test2PhantomSD::ProcessHits(G4Step * aStep, G4TouchableHistory *) {

  // total energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double trklen = aStep->GetStepLength();
  if(verboseLevel > 1) G4cout << "Next step edep [MeV] = " << edep/MeV << G4endl;
  if(edep > 0. || trklen > 0.) {

    G4TouchableHistory * hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* physVol = hist->GetVolume();
    G4int copyIDinX = hist->GetReplicaNumber(0);
    G4int copyIDinY = hist->GetReplicaNumber(1);
    G4int copyIDinZ = hist->GetReplicaNumber(2);

    if(fCellID[copyIDinX][copyIDinY][copyIDinZ] == -1) {
      Test2PhantomHit* phantomHit
	= new Test2PhantomHit(physVol->GetLogicalVolume(), copyIDinX, copyIDinY, copyIDinZ);
      phantomHit->SetEdep(edep);
      phantomHit->SetTrackLength(trklen);
      phantomHit->SetParticleName(aStep->GetTrack()->GetParticleDefinition()->GetParticleName());
      G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
      aTrans.Invert();
      phantomHit->SetPos(aTrans.NetTranslation());
      phantomHit->SetRot(aTrans.NetRotation());
      G4int icell = fPhantomCollection->insert(phantomHit);
      fCellID[copyIDinX][copyIDinY][copyIDinZ] = icell - 1;
      if(verboseLevel > 0) {
	G4cout << " New Calorimeter Hit on fCellID " 
	       << copyIDinX << ", " << copyIDinY << ", " << copyIDinZ << G4endl; }
    } else { 
      (*fPhantomCollection)[fCellID[copyIDinX][copyIDinY][copyIDinZ]]->AddEdep(edep);
      (*fPhantomCollection)[fCellID[copyIDinX][copyIDinY][copyIDinZ]]->AddTrackLength(trklen);
      if(verboseLevel > 0) {
	G4cout << " Energy added to fCellID " 
	       << copyIDinX << ", " << copyIDinY << ", " << copyIDinZ << G4endl; }
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

