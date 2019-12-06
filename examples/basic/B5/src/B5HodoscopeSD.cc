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
/// \file B5HodoscopeSD.cc
/// \brief Implementation of the B5HodoscopeSD class

#include "B5HodoscopeSD.hh"
#include "B5HodoscopeHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5HodoscopeSD::B5HodoscopeSD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert( "hodoscopeColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5HodoscopeSD::~B5HodoscopeSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5HodoscopeSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new B5HodoscopeHitsCollection
  (SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) { 
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool B5HodoscopeSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  
  auto preStepPoint = step->GetPreStepPoint();
  auto touchable = preStepPoint->GetTouchable();
  auto copyNo = touchable->GetVolume()->GetCopyNo();
  auto hitTime = preStepPoint->GetGlobalTime();
  
  // check if this finger already has a hit
  auto ix = -1;
  for (std::size_t i=0;i<fHitsCollection->entries();++i) {
    if ((*fHitsCollection)[i]->GetID()==copyNo) {
      ix = i;
      break;
    }
  }

  if (ix>=0) {
    // if it has, then take the earlier time
    if ((*fHitsCollection)[ix]->GetTime()>hitTime) { 
      (*fHitsCollection)[ix]->SetTime(hitTime); 
    }
  }
  else {
    // if not, create a new hit and set it to the collection
    auto hit = new B5HodoscopeHit(copyNo,hitTime);
    auto physical = touchable->GetVolume();
    hit->SetLogV(physical->GetLogicalVolume());
    auto transform = touchable->GetHistory()->GetTopTransform();
    transform.Invert();
    hit->SetRot(transform.NetRotation());
    hit->SetPos(transform.NetTranslation());
    fHitsCollection->insert(hit);
  }    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
