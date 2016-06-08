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

#include "ExN04CalorimeterSD.hh"
#include "ExN04CalorimeterHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"

ExN04CalorimeterSD::ExN04CalorimeterSD(G4String name)
:G4VSensitiveDetector(name),
 numberOfCellsInZ(20),numberOfCellsInPhi(48)
{
  G4String HCname;
  collectionName.insert(HCname="calCollection");
}

ExN04CalorimeterSD::~ExN04CalorimeterSD()
{;}

void ExN04CalorimeterSD::Initialize(G4HCofThisEvent*HCE)
{
  CalCollection = new ExN04CalorimeterHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for(int j=0;j<numberOfCellsInZ;j++)
  for(int k=0;k<numberOfCellsInPhi;k++)
  {
    CellID[j][k] = -1;
  }
}

G4bool ExN04CalorimeterSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  if(!ROhist) return false;
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;

  G4VPhysicalVolume* physVol = ROhist->GetVolume();
  //ROhist->MoveUpHistory();
  G4VPhysicalVolume* mothVol = ROhist->GetVolume(1);
  int copyIDinZ = ROhist->GetReplicaNumber();
  int copyIDinPhi = ROhist->GetReplicaNumber(1);

  if(CellID[copyIDinZ][copyIDinPhi]==-1)
  {
    ExN04CalorimeterHit* calHit
      = new ExN04CalorimeterHit
            (physVol->GetLogicalVolume(),copyIDinZ,copyIDinPhi);
    G4RotationMatrix rotM;
    if(physVol->GetObjectRotation()) rotM = *(physVol->GetObjectRotation());
    calHit->SetEdep( edep );
    calHit->SetPos( physVol->GetTranslation() );
    calHit->SetRot( rotM );
    int icell = CalCollection->insert( calHit );
    CellID[copyIDinZ][copyIDinPhi] = icell - 1;
    if(verboseLevel>0)
    { G4cout << " New Calorimeter Hit on CellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }
  else
  { 
    (*CalCollection)[CellID[copyIDinZ][copyIDinPhi]]->AddEdep(edep);
    if(verboseLevel>0)
    { G4cout << " Energy added to CellID " 
           << copyIDinZ << " " << copyIDinPhi << G4endl; }
  }

  return true;
}

void ExN04CalorimeterSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, CalCollection );
}

void ExN04CalorimeterSD::clear()
{
} 

void ExN04CalorimeterSD::DrawAll()
{
} 

void ExN04CalorimeterSD::PrintAll()
{
} 

