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
// $Id: Tst23SensitiveDetector.cc,v 1.1 2001-12-14 14:53:43 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst23SensitiveDetector.hh"
#include "Tst23Hit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

Tst23SensitiveDetector::Tst23SensitiveDetector(G4String name)
  :G4VSensitiveDetector(name)
{
  G4String HCname="HitsCollection";
  collectionName.insert(HCname);
}

Tst23SensitiveDetector::~Tst23SensitiveDetector()
{
}

void Tst23SensitiveDetector::Initialize(G4HCofThisEvent*HCE)
{
  HitCollection = new Tst23HitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  for (size_t idx=0; idx<2; idx++) {
    Tst23Hit* hit = new Tst23Hit(idx);
    hit->SetEdep( 0.0 );
    HitCollection->insert( hit );
  }
}

G4bool Tst23SensitiveDetector::ProcessHits(G4Step*aStep,
					    G4TouchableHistory*ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep<=0.) return false;

  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();
  int copyID = physVol->GetCopyNo();
  
  (*HitCollection)[copyID]->AddEdep( edep );

  if(verboseLevel>0) {
    G4cout << " Energy added to CellID " << copyID << G4endl; 
  }

  return true;
}

void Tst23SensitiveDetector::EndOfEvent(G4HCofThisEvent*HCE)
{
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, HitCollection );
}

void Tst23SensitiveDetector::clear()
{
} 

void Tst23SensitiveDetector::DrawAll()
{
} 

void Tst23SensitiveDetector::PrintAll()
{
} 







