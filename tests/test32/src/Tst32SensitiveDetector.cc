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
// $Id: Tst32SensitiveDetector.cc,v 1.1 2002-06-13 12:16:36 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst32SensitiveDetector.hh"
#include "Tst32Hit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

Tst32SensitiveDetector::Tst32SensitiveDetector(G4String name)
  :G4VSensitiveDetector(name)
{
  G4String HCname="HitsCollection";
  collectionName.insert(HCname);
}

Tst32SensitiveDetector::~Tst32SensitiveDetector()
{
}

void Tst32SensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
  HitCollection = new Tst32HitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  for (size_t idx=0; idx<10; idx++) {
    Tst32Hit* hit = new Tst32Hit(idx);
    hit->SetEdep( 0.0 );
    HitCollection->insert( hit );
  }
}

G4bool Tst32SensitiveDetector::ProcessHits(G4Step*aStep,
					    G4TouchableHistory*ROhist)
{

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4Track * aTrack = aStep->GetTrack();
  G4double weight = aTrack->GetWeight();
  if(weight < DBL_MIN || weight > DBL_MAX || edep < 0) return false;

  edep *= weight;

  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int copyID = physVol->GetCopyNo();
  
  (*HitCollection)[copyID]->AddEdep( edep );

  if(verboseLevel>0) {
    G4cout << " Energy added to CellID " << copyID << G4endl; 
  }

  return true;
}

void Tst32SensitiveDetector::EndOfEvent(G4HCofThisEvent*HCE)
{
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, HitCollection );
}

void Tst32SensitiveDetector::clear()
{
} 

void Tst32SensitiveDetector::DrawAll()
{
} 

void Tst32SensitiveDetector::PrintAll()
{
} 
