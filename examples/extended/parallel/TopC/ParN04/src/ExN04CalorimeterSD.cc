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
/// \file ExN04CalorimeterSD.cc
/// \brief Implementation of the ExN04CalorimeterSD class
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
#include "G4SystemOfUnits.hh"
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

void ExN04CalorimeterSD::Initialize(G4HCofThisEvent*)
{
  CalCollection = new ExN04CalorimeterHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  for(G4int j=0;j<numberOfCellsInZ;j++)
  for(G4int k=0;k<numberOfCellsInPhi;k++)
  {
    CellID[j][k] = -1;
  }
  verboseLevel = 0;
}

G4bool ExN04CalorimeterSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
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

  if(CellID[copyIDinZ][copyIDinPhi]==-1)
  {
    ExN04CalorimeterHit* calHit
      = new ExN04CalorimeterHit
            (physVol->GetLogicalVolume(),copyIDinZ,copyIDinPhi);
    calHit->SetEdep( edep );
    G4AffineTransform aTrans = ROhist->GetHistory()->GetTopTransform();
    aTrans.Invert();
    calHit->SetPos(aTrans.NetTranslation());
    calHit->SetRot(aTrans.NetRotation());
    G4int icell = CalCollection->insert( calHit );
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

