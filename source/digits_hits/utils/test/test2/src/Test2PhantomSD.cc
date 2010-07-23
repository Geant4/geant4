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
   numberOfCellsInX(10), numberOfCellsInY(10), numberOfCellsInZ(10)
 {

  G4String HCname;
  collectionName.insert(HCname = "phantomCollection");

}

Test2PhantomSD::~Test2PhantomSD() {
  ;
}

void Test2PhantomSD::Initialize(G4HCofThisEvent *) {

  phantomCollection = new Test2PhantomHitsCollection(SensitiveDetectorName,
						    collectionName[0]); 
  for(G4int i = 0; i < numberOfCellsInZ; i++) {
    for(G4int j = 0; j < numberOfCellsInY; j++) {
      for(G4int k = 0; k < numberOfCellsInX; k++) {
	CellID[i][j][k] = -1;
      }
    }
  }
  verboseLevel = 0;
}

G4bool Test2PhantomSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  if(!ROhist) return false;
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(verboseLevel>1) G4cout << "Next step edep(MeV) = " << edep/MeV << G4endl;
  if(edep==0.) return false;

  G4VPhysicalVolume* physVol = ROhist->GetVolume();
  //ROhist->MoveUpHistory();
  //G4VPhysicalVolume* mothVol = ROhist->GetVolume(1);
  G4int copyIDinZ = ROhist->GetReplicaNumber();
  G4int copyIDinPhi = ROhist->GetReplicaNumber(1);

  if(CellID[copyIDinZ][copyIDinPhi][0]==-1)
  {
    Test2PhantomHit* calHit
      = new Test2PhantomHit(physVol->GetLogicalVolume(),copyIDinZ,copyIDinPhi);
    calHit->SetEdep( edep );
    G4AffineTransform aTrans = ROhist->GetHistory()->GetTopTransform();
    aTrans.Invert();
    calHit->SetPos(aTrans.NetTranslation());
    calHit->SetRot(aTrans.NetRotation());
    G4int icell = phantomCollection->insert( calHit );
    CellID[copyIDinZ][copyIDinPhi][0] = icell - 1;
    if(verboseLevel>0)
    { G4cout << " New Calorimeter Hit on CellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }
  else
  { 
    (*phantomCollection)[CellID[copyIDinZ][copyIDinPhi][0]]->AddEdep(edep);
    if(verboseLevel>0)
    { G4cout << " Energy added to CellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }

  return true;
}

void Test2PhantomSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, phantomCollection );
}

void Test2PhantomSD::clear()
{
} 

void Test2PhantomSD::DrawAll()
{
} 

void Test2PhantomSD::PrintAll()
{
} 

