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
/// \file optical/wls/src/WLSPhotonDetSD.cc
/// \brief Implementation of the WLSPhotonDetSD class
//
//
#include "WLSPhotonDetSD.hh"

#include "WLSPhotonDetHit.hh"
#include "WLSUserTrackInformation.hh"

#include "G4OpticalPhoton.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhotonDetSD::WLSPhotonDetSD(G4String name)
  : G4VSensitiveDetector(name)
  , fPhotonDetHitCollection(nullptr)
  , fHCID(-1)
{
  collectionName.insert("PhotonDetHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WLSPhotonDetSD::~WLSPhotonDetSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WLSPhotonDetSD::Initialize(G4HCofThisEvent* HCE)
{
  fPhotonDetHitCollection =
    new WLSPhotonDetHitsCollection(SensitiveDetectorName, collectionName[0]);

  if(fHCID < 0)
  {
    fHCID =
      G4SDManager::GetSDMpointer()->GetCollectionID(fPhotonDetHitCollection);
  }
  HCE->AddHitsCollection(fHCID, fPhotonDetHitCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool WLSPhotonDetSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  if(!aStep)
    return false;
  G4Track* theTrack = aStep->GetTrack();

  // Need to know if this is an optical photon
  if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
  {
    return false;
  }

  // Find out information regarding the hit
  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  WLSUserTrackInformation* trackInformation =
    (WLSUserTrackInformation*) theTrack->GetUserInformation();

  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*) (thePostPoint->GetTouchable());

  G4ThreeVector photonExit   = trackInformation->GetExitPosition();
  G4ThreeVector photonArrive = thePostPoint->GetPosition();
  G4double arrivalTime       = theTrack->GetGlobalTime();
  G4double energy            = theTrack->GetTotalEnergy();

  // Convert the global coordinate for arriving photons into
  // the local coordinate of the detector
  photonArrive =
    theTouchable->GetHistory()->GetTopTransform().TransformPoint(photonArrive);

  // Creating the hit and add it to the collection
  fPhotonDetHitCollection->insert(
    new WLSPhotonDetHit(photonExit, photonArrive, arrivalTime, energy));

  return true;
}

void WLSPhotonDetSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) {
     G4int nofHits = fPhotonDetHitCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event there are " << nofHits
            << " hits in the photon detector: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fPhotonDetHitCollection)[i]->Print();
  }
}

