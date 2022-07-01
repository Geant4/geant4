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
#include "Par04SensitiveDetector.hh"
#include <CLHEP/Vector/Rotation.h>     // for HepRotation
#include <CLHEP/Vector/ThreeVector.h>  // for Hep3Vector
#include <cmath>                       // for floor
#include <G4CollectionNameVector.hh>   // for G4CollectionNameVector
#include <G4FastHit.hh>                // for G4FastHit
#include <G4FastTrack.hh>              // for G4FastTrack
#include <G4RotationMatrix.hh>         // for G4RotationMatrix
#include <G4StepPoint.hh>              // for G4StepPoint
#include <G4THitsCollection.hh>        // for G4THitsCollection
#include <G4ThreeVector.hh>            // for G4ThreeVector
#include <G4VSensitiveDetector.hh>     // for G4VSensitiveDetector
#include <G4VUserEventInformation.hh>  // for G4VUserEventInformation
#include <cstddef>                     // for size_t
#include <vector>                      // for vector
#include "G4Event.hh"                  // for G4Event
#include "G4EventManager.hh"           // for G4EventManager
#include "G4HCofThisEvent.hh"          // for G4HCofThisEvent
#include "G4SDManager.hh"              // for G4SDManager
#include "G4Step.hh"                   // for G4Step
#include "G4Track.hh"                  // for G4Track
#include "Par04EventInformation.hh"    // for Par04EventInformation
#include "Par04Hit.hh"                 // for Par04Hit, Par04HitsCollection

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04SensitiveDetector::Par04SensitiveDetector(G4String aName)
  : G4VSensitiveDetector(aName)
{
  collectionName.insert("hits");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04SensitiveDetector::Par04SensitiveDetector(G4String aName, G4ThreeVector aNb,
                                               G4ThreeVector aSize)
  : G4VSensitiveDetector(aName)
  , fMeshNbOfCells(aNb)
  , fMeshSizeOfCells(aSize)
{
  collectionName.insert("hits");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04SensitiveDetector::~Par04SensitiveDetector() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Par04SensitiveDetector::Initialize(G4HCofThisEvent* aHCE)
{
  fHitsCollection = new Par04HitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHitCollectionID < 0)
  {
    fHitCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  aHCE->AddHitsCollection(fHitCollectionID, fHitsCollection);

  // fill calorimeter hits with zero energy deposition
  for(G4int iphi = 0; iphi < fMeshNbOfCells.y(); iphi++)
    for(G4int irho = 0; irho < fMeshNbOfCells.x(); irho++)
      for(G4int iz = 0; iz < fMeshNbOfCells.z(); iz++)
      {
        Par04Hit* hit = new Par04Hit();
        fHitsCollection->insert(hit);
      }
  // reset entrance position
  fEntrancePosition.set(-1, -1, -1);
  fEntranceDirection.set(-1, -1, -1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.)
    return true;

  auto hit = RetrieveAndSetupHit(aStep->GetPostStepPoint()->GetPosition());
  if(hit == nullptr)
    return true;

  // Add energy deposit from G4Step
  hit->AddEdep(edep);

  // Fill time information from G4Step
  // If it's already filled, choose hit with earliest global time
  if(hit->GetTime() == -1 || hit->GetTime() > aStep->GetTrack()->GetGlobalTime())
    hit->SetTime(aStep->GetTrack()->GetGlobalTime());

  // Set hit type to full simulation (only if hit is not already marked as fast
  // sim)
  if(hit->GetType() != 1)
    hit->SetType(0);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Par04SensitiveDetector::ProcessHits(const G4FastHit* aHit, const G4FastTrack* aTrack,
                                           G4TouchableHistory*)
{
  G4double edep = aHit->GetEnergy();
  if(edep == 0.)
    return true;

  auto hit = RetrieveAndSetupHit(aHit->GetPosition());
  if(hit == nullptr)
    return true;

  // Add energy deposit from G4FastHit
  hit->AddEdep(edep);

  // Fill time information from G4FastTrack
  // If it's already filled, choose hit with earliest global time
  if(hit->GetTime() == -1 || hit->GetTime() > aTrack->GetPrimaryTrack()->GetGlobalTime())
  {
    hit->SetTime(aTrack->GetPrimaryTrack()->GetGlobalTime());
  }

  // Set hit type to fast simulation (even if hit was already marked as full
  // sim, overwrite it)
  hit->SetType(1);

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par04Hit* Par04SensitiveDetector::RetrieveAndSetupHit(G4ThreeVector aGlobalPosition)
{
  if(fEntrancePosition.x() == -1)
  {
    Par04EventInformation* info = dynamic_cast<Par04EventInformation*>(
      G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation());
    if(info == nullptr)
      return nullptr;
    fEntrancePosition  = info->GetPosition();
    fEntranceDirection = info->GetDirection();
  }

  auto delta = aGlobalPosition - fEntrancePosition;

  // Calculate rotation matrix along the particle momentum direction
  // It will rotate the shower axes to match the incoming particle direction
  G4RotationMatrix rotMatrix = G4RotationMatrix();
  double particleTheta       = fEntranceDirection.theta();
  double particlePhi         = fEntranceDirection.phi();
  rotMatrix.rotateZ(-particlePhi);
  rotMatrix.rotateY(-particleTheta);
  G4RotationMatrix rotMatrixInv = CLHEP::inverseOf(rotMatrix);

  delta = rotMatrix * delta;

  G4int rhoNo = std::floor(delta.perp() / fMeshSizeOfCells.x());
  G4int phiNo = std::floor((CLHEP::pi + delta.phi()) / fMeshSizeOfCells.y());
  G4int zNo   = std::floor(delta.z() / fMeshSizeOfCells.z());

  std::size_t hitID =
    fMeshNbOfCells.x() * fMeshNbOfCells.z() * phiNo + fMeshNbOfCells.z() * rhoNo + zNo;

  if(hitID >= fHitsCollection->entries() || zNo >= fMeshNbOfCells.z() ||
     rhoNo >= fMeshNbOfCells.x() || zNo < 0)
  {
    return nullptr;
  }

  Par04Hit* hit = (*fHitsCollection)[hitID];

  if(hit->GetRhoId() < 0)
  {
    hit->SetRhoId(rhoNo);
    hit->SetPhiId(phiNo);
    hit->SetZid(zNo);
    hit->SetRot(rotMatrixInv);
    hit->SetPos(fEntrancePosition +
                rotMatrixInv * G4ThreeVector(0, 0, (zNo + 0.5) * fMeshSizeOfCells.z()));
  }
  return hit;
}
