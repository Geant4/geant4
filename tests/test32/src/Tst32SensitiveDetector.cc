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
// $Id: Tst32SensitiveDetector.cc,v 1.5 2006-06-29 21:59:23 gunter Exp $
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
  : G4VSensitiveDetector(name), HitCollection(0)
{
  G4String HCname="HitsCollection";
  collectionName.insert(HCname);
}

Tst32SensitiveDetector::~Tst32SensitiveDetector()
{
}

void Tst32SensitiveDetector::Initialize(G4HCofThisEvent*)
{
  HitCollection = new Tst32HitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  for (size_t idx=0; idx<10; idx++)
  {
    Tst32Hit* hit = new Tst32Hit(idx);
    hit->SetEdep( 0.0 );
    HitCollection->insert( hit );
  }
}

G4bool Tst32SensitiveDetector::ProcessHits(G4Step* aStep,
					   G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4Track * aTrack = aStep->GetTrack();
  G4double weight = aTrack->GetWeight();
  if(weight < DBL_MIN || weight > DBL_MAX || edep < 0) return false;

  edep *= weight;

  G4int copyID = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber();
  (*HitCollection)[copyID%10]->AddEdep( edep );

  if (verboseLevel>0)
    G4cout << " Energy of: " << edep
           << " MeV, added to CellID " << copyID%10 << G4endl;

  return true;
}

void Tst32SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE)
{
  static G4int HCID = -1;
  if (HCID<0)
  {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
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
