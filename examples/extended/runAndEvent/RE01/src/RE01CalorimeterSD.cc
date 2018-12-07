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
/// \file runAndEvent/RE01/src/RE01CalorimeterSD.cc
/// \brief Implementation of the RE01CalorimeterSD class
//
//

#include "RE01CalorimeterSD.hh"
#include "RE01CalorimeterHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterSD::RE01CalorimeterSD(G4String name)
  :G4VSensitiveDetector(name), fCalCollection(0),
   fNumberOfCellsInZ(20),fNumberOfCellsInPhi(48)
{
  // Initialize data member.
  for(G4int j=0;j<fNumberOfCellsInZ;j++)
  for(G4int k=0;k<fNumberOfCellsInPhi;k++)
  {
    fCellID[j][k] = -1;
  }
  //
  G4String HCname;
  collectionName.insert(HCname="calCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE01CalorimeterSD::~RE01CalorimeterSD()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE01CalorimeterSD::Initialize(G4HCofThisEvent* HCE)
{
  fCalCollection = new RE01CalorimeterHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for(G4int j=0;j<fNumberOfCellsInZ;j++)
  for(G4int k=0;k<fNumberOfCellsInPhi;k++)
  {
    fCellID[j][k] = -1;
  }

  static G4int HCID = -1;
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, fCalCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool RE01CalorimeterSD::ProcessHits(G4Step*aStep,G4TouchableHistory*)
{
//***** RE05CalorimeterSD has been migrated to Geant4 version 10 that does not
//***** support Readout Geometry in multi-threaded mode. Now RE05CalorimeterSD
//***** is assigned to a dedicaed parallel world. The pointer "aStep" points to
//***** a G4Step object for the parallel world.

  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

  const G4VTouchable* ROhist = aStep->GetPreStepPoint()->GetTouchable();
  G4int copyIDinZ = ROhist->GetReplicaNumber();
  G4int copyIDinPhi = ROhist->GetReplicaNumber(1);

  if(fCellID[copyIDinZ][copyIDinPhi]==-1)
  {
    RE01CalorimeterHit* calHit = new RE01CalorimeterHit
            (ROhist->GetVolume()->GetLogicalVolume(),copyIDinZ,copyIDinPhi);
    calHit->SetEdep( edep );
    G4AffineTransform aTrans = ROhist->GetHistory()->GetTopTransform();
    aTrans.Invert();
    calHit->SetPos(aTrans.NetTranslation());
    calHit->SetRot(aTrans.NetRotation());
    calHit->SetTrackInformation(aStep->GetTrack());
    G4int icell = fCalCollection->insert( calHit );
    fCellID[copyIDinZ][copyIDinPhi] = icell - 1;
    if(verboseLevel>0)
    { G4cout << " New Calorimeter Hit on CellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }
  else
  { 
    (*fCalCollection)[fCellID[copyIDinZ][copyIDinPhi]]->AddEdep(edep);
    (*fCalCollection)[fCellID[copyIDinZ][copyIDinPhi]]
                             ->SetTrackInformation(aStep->GetTrack());
    if(verboseLevel>0)
    { G4cout << " Energy added to CellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }

  return true;
}

