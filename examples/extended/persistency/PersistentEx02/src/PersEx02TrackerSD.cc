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
// $Id: PersEx02TrackerSD.cc,v 1.6 2001/07/11 09:58:16 gunter Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//

#include <assert.h>
#include "PersEx02TrackerSD.hh"
#include "PersEx02TrackerHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4PersistentHitMan.hh"
#include "G4PHCofThisEvent.hh"
#include "G4ios.hh"

PersEx02TrackerSD::PersEx02TrackerSD(G4String name)
 : G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="EvenCollection");
  collectionName.insert(HCname="OddCollection");

  f_hitMan = G4PersistentHitMan::GetPersistentHitMan();
  hitContainer = f_hitMan->GetContainer();
  assert(hitContainer != 0);
}

PersEx02TrackerSD::~PersEx02TrackerSD(){;}

void PersEx02TrackerSD::Initialize(G4HCofThisEvent*HCE)
{
  HepRef(G4PHCofThisEvent) aPHC = new( hitContainer ) G4PHCofThisEvent;
  f_hitMan->SetCurrentPHCofThisEvent(aPHC);

  EvenCollection = new( hitContainer )
        PersEx02TrackerHitsCollection(SensitiveDetectorName,collectionName[0]); 
  OddCollection  = new( hitContainer )
        PersEx02TrackerHitsCollection(SensitiveDetectorName,collectionName[1]); 
}

G4bool PersEx02TrackerSD::ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if( edep == 0. ) return true;

  const G4VPhysicalVolume* physVol 
    = aStep->GetPreStepPoint()->GetPhysicalVolume();

  HepRef(PersEx02TrackerHit) newHit = new( hitContainer ) PersEx02TrackerHit();
  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() );

  if( physVol->GetName() == "TPhys1" )
  { EvenCollection->Insert( newHit ); }
  else if( physVol->GetName() == "TPhys2" )
  { OddCollection->Insert( newHit ); }

  return true;
}

void PersEx02TrackerSD::EndOfEvent(G4HCofThisEvent*HCE)
{
  static G4int IDe = -1;
  static G4int IDo = -1;
  if(IDe<0)
  { IDe = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  if(IDo<0)
  { IDo = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[1]); }

  // Shrink the array size to the actual hit entries
  EvenCollection->ShrinkToFit();
  OddCollection->ShrinkToFit();

  // Register this HitCollection to Persistent HCofThisEvent
  G4PHCofThisEvent* aPHC = f_hitMan->GetCurrentPHCofThisEvent();
  aPHC->AddHitsCollection( IDe, EvenCollection ); 
  aPHC->AddHitsCollection( IDo,  OddCollection ); 
}

void PersEx02TrackerSD::clear()
{
} 

void PersEx02TrackerSD::DrawAll()
{
} 

void PersEx02TrackerSD::PrintAll()
{
} 

