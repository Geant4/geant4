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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SpecializedTrackingManager.hh"

#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4RegionStore.hh"
#include "G4StackManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrackingManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SpecializedTrackingManager::SpecializedTrackingManager() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SpecializedTrackingManager::~SpecializedTrackingManager() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SpecializedTrackingManager::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if (fBackRegion == nullptr) {
    fBackRegion = G4RegionStore::GetInstance()->GetRegion("Back", false);
  }

  G4ProcessManager* pManager = part.GetProcessManager();
  G4ProcessManager* pManagerShadow = part.GetMasterProcessManager();

  G4ProcessVector* pVector = pManager->GetProcessList();
  for (std::size_t j = 0; j < pVector->size(); ++j) {
    if (pManagerShadow == pManager) {
      (*pVector)[j]->BuildPhysicsTable(part);
    }
    else {
      (*pVector)[j]->BuildWorkerPhysicsTable(part);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SpecializedTrackingManager::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  G4ProcessManager* pManager = part.GetProcessManager();
  G4ProcessManager* pManagerShadow = part.GetMasterProcessManager();

  G4ProcessVector* pVector = pManager->GetProcessList();
  for (std::size_t j = 0; j < pVector->size(); ++j) {
    if (pManagerShadow == pManager) {
      (*pVector)[j]->PreparePhysicsTable(part);
    }
    else {
      (*pVector)[j]->PrepareWorkerPhysicsTable(part);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SpecializedTrackingManager::HandOverOneTrack(G4Track* aTrack)
{
  if (aTrack->GetKineticEnergy() < 100 * MeV) {
    // If the particle energy is lower than 100 MeV, track it immediately by
    // passing to the generic G4TrackingManager. This avoids storing lower
    // energy particles in the buffer and feeding it through the specialized
    // tracking.
    G4EventManager* eventManager = G4EventManager::GetEventManager();
    G4TrackingManager* trackManager = eventManager->GetTrackingManager();

    trackManager->ProcessOneTrack(aTrack);
    if (aTrack->GetTrackStatus() != fStopAndKill) {
      G4Exception("SpecializedTrackingManager::HandOverOneTrack", "NotStopped", FatalException,
        "track was not stopped");
    }

    G4TrackVector* secondaries = trackManager->GimmeSecondaries();
    eventManager->StackTracks(secondaries);
    delete aTrack;
    return;
  }

  fBufferedTracks.push_back(aTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SpecializedTrackingManager::FlushEvent()
{
  G4EventManager* eventManager = G4EventManager::GetEventManager();
  G4TrackingManager* trackManager = eventManager->GetTrackingManager();
  G4SteppingManager* steppingManager = trackManager->GetSteppingManager();
  G4TrackVector* secondaries = trackManager->GimmeSecondaries();

  for (G4Track* aTrack : fBufferedTracks) {
    // Clear secondary particle vector
    for (std::size_t itr = 0; itr < secondaries->size(); ++itr) {
      delete (*secondaries)[itr];
    }
    secondaries->clear();

    steppingManager->SetInitialStep(aTrack);

    G4UserTrackingAction* userTrackingAction = trackManager->GetUserTrackingAction();
    if (userTrackingAction != nullptr) {
      userTrackingAction->PreUserTrackingAction(aTrack);
    }

    // Give SteppingManger the maxmimum number of processes
    steppingManager->GetProcessNumber();

    // Give track the pointer to the Step
    aTrack->SetStep(steppingManager->GetStep());

    // Inform beginning of tracking to physics processes
    aTrack->GetDefinition()->GetProcessManager()->StartTracking(aTrack);

    // Track the particle Step-by-Step while it is alive
    while ((aTrack->GetTrackStatus() == fAlive) || (aTrack->GetTrackStatus() == fStopButAlive)) {
      G4Region* region = aTrack->GetVolume()->GetLogicalVolume()->GetRegion();
      if (region == fBackRegion) {
        StepInBackRegion(aTrack);
      }
      else {
        StepOutside(aTrack);
      }
    }

    aTrack->GetDefinition()->GetProcessManager()->EndTracking();

    if (userTrackingAction != nullptr) {
      userTrackingAction->PostUserTrackingAction(aTrack);
    }

    eventManager->StackTracks(secondaries);
    delete aTrack;
  }

  fBufferedTracks.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SpecializedTrackingManager::StepInBackRegion(G4Track* aTrack)
{
  G4EventManager* eventManager = G4EventManager::GetEventManager();
  G4TrackingManager* trackManager = eventManager->GetTrackingManager();
  G4SteppingManager* steppingManager = trackManager->GetSteppingManager();

  // Track the particle Step-by-Step while it is alive and inside the "Back"
  // region of the detector. Implement a low-energy cut-off for particles
  // below 100 MeV. More specialized handling would also be possible, such
  // as only killing particles in non-sensitive materials / volumes.
  while ((aTrack->GetTrackStatus() == fAlive) || (aTrack->GetTrackStatus() == fStopButAlive)) {
    aTrack->IncrementCurrentStepNumber();
    steppingManager->Stepping();

    if (aTrack->GetTrackStatus() != fStopAndKill) {
      // Switch the touchable to update the volume, which is checked in the
      // condition below and at the call site.
      aTrack->SetTouchableHandle(aTrack->GetNextTouchableHandle());
      G4Region* region = aTrack->GetVolume()->GetLogicalVolume()->GetRegion();
      if (region != fBackRegion) {
        return;
      }

      if (aTrack->GetKineticEnergy() < 100 * MeV) {
        // Kill the particle.
        aTrack->SetTrackStatus(fStopAndKill);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SpecializedTrackingManager::StepOutside(G4Track* aTrack)
{
  G4EventManager* eventManager = G4EventManager::GetEventManager();
  G4TrackingManager* trackManager = eventManager->GetTrackingManager();
  G4SteppingManager* steppingManager = trackManager->GetSteppingManager();

  // Track the particle Step-by-Step while it is alive and still outside of
  // the "Back" region.
  while ((aTrack->GetTrackStatus() == fAlive) || (aTrack->GetTrackStatus() == fStopButAlive)) {
    aTrack->IncrementCurrentStepNumber();
    steppingManager->Stepping();

    if (aTrack->GetTrackStatus() != fStopAndKill) {
      // Switch the touchable to update the volume, which is checked in the
      // condition below and at the call site.
      aTrack->SetTouchableHandle(aTrack->GetNextTouchableHandle());
      G4Region* region = aTrack->GetVolume()->GetLogicalVolume()->GetRegion();
      if (region == fBackRegion) {
        return;
      }
    }
  }
}
