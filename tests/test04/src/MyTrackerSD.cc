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

#include "MyTrackerSD.hh"
#include "MyTrackerHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

MyTrackerSD::MyTrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="EvenCollection");
  collectionName.insert(HCname="OddCollection");
}

MyTrackerSD::~MyTrackerSD(){;}

void MyTrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  EvenCollection = new 
        MyTrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
  OddCollection = new 
        MyTrackerHitsCollection(SensitiveDetectorName,collectionName[1]); 
}

G4bool MyTrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();

  MyTrackerHit* newHit = new MyTrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );
  if( physVol->GetName() == "TPhys1" )
  { EvenCollection->insert( newHit ); }
  else if( physVol->GetName() == "TPhys2" )
  { OddCollection->insert( newHit ); }
  return true;
}

void MyTrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int IDe = -1;
  static G4int IDo = -1;
  if(IDe<0)
  { IDe = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  if(IDo<0)
  { IDo = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]); }
  HCE->AddHitsCollection( IDe, EvenCollection );
  HCE->AddHitsCollection( IDo, OddCollection ); 
}

void MyTrackerSD::clear()
{
} 

void MyTrackerSD::DrawAll()
{
} 

void MyTrackerSD::PrintAll()
{
} 

