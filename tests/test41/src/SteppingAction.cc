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
#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* adet, RunAction* arun)
:G4UserSteppingAction(),det(adet),run(arun) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4Track* track = aStep->GetTrack();
  if(track->GetTrackID() != 1) {
    track->SetTrackStatus(fStopAndKill);
    return;
  }
  G4VPhysicalVolume* vol1 = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* vol2 = aStep->GetPostStepPoint()->GetPhysicalVolume();
  if(vol1 == det->GetPhysAbsorber() && vol2 == det->GetPhysWorld()) {
    G4ThreeVector dir = track->GetMomentumDirection();
    /*
    G4cout << "Next track vx= " << std::abs(dir.x())
	   << " vy= " <<  std::abs(dir.y()) << G4endl;
    */
    if(dir.z() > 0.0) {
      run->AddProjectileTheta(std::asin(std::abs(dir.y())));
      track->SetTrackStatus(fStopAndKill);
    }
  }  
}

