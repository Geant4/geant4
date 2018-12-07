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
/// \file RE05/src/RE05CalorimeterSD.cc
/// \brief Implementation of the RE05CalorimeterSD class
//

#include "RE05CalorimeterSD.hh"
#include "RE05CalorimeterHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05CalorimeterSD::RE05CalorimeterSD(G4String name)
:G4VSensitiveDetector(name),
 fNumberOfCellsInZ(20),fNumberOfCellsInPhi(48)
{
  G4String HCname;
  collectionName.insert(HCname="calCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RE05CalorimeterSD::~RE05CalorimeterSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05CalorimeterSD::Initialize(G4HCofThisEvent*)
{
  fCalCollection = new RE05CalorimeterHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for(G4int j=0;j<fNumberOfCellsInZ;j++)
  for(G4int k=0;k<fNumberOfCellsInPhi;k++)
  {
    fCellID[j][k] = -1;
  }
  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool RE05CalorimeterSD::ProcessHits(G4Step*aStep,G4TouchableHistory*)
{
//***** RE05CalorimeterSD has been migrated to Geant4 version 10 that does not
//***** support Readout Geometry in multi-threaded mode. Now RE05CalorimeterSD
//***** is assigned to a dedicated parallel world. The pointer "aStep" points to
//***** a G4Step object for the parallel world.

  G4double edep = aStep->GetTotalEnergyDeposit();
  if(verboseLevel>1) G4cout << "Next step edep(MeV) = " << edep/MeV << G4endl;
  if(edep==0.) return false;

  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int copyIDinZ = touchable->GetReplicaNumber(0);
  G4int copyIDinPhi = touchable->GetReplicaNumber(1);

  if(fCellID[copyIDinZ][copyIDinPhi]==-1)
  {
    RE05CalorimeterHit* calHit
      = new RE05CalorimeterHit
          (touchable->GetVolume()->GetLogicalVolume(),copyIDinZ,copyIDinPhi);
    calHit->SetEdep( edep );
    G4AffineTransform aTrans = touchable->GetHistory()->GetTopTransform();
    aTrans.Invert();
    calHit->SetPos(aTrans.NetTranslation());
    calHit->SetRot(aTrans.NetRotation());
    G4int icell = fCalCollection->insert( calHit );
    fCellID[copyIDinZ][copyIDinPhi] = icell - 1;
    if(verboseLevel>0)
    { G4cout << " New Calorimeter Hit on fCellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }
  else
  { 
    (*fCalCollection)[fCellID[copyIDinZ][copyIDinPhi]]->AddEdep(edep);
    if(verboseLevel>0)
    { G4cout << " Energy added to fCellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05CalorimeterSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, fCalCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05CalorimeterSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05CalorimeterSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RE05CalorimeterSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
