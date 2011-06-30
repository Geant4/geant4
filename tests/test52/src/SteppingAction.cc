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
#include "RegionInformation.hh"
#include "TrackInformation.hh"
#include "AnalysisManager.hh"
#include "G4Step.hh"


SteppingAction::SteppingAction() {

}


SteppingAction::~SteppingAction() {

}


void SteppingAction::UserSteppingAction(const G4Step* step) {

  if(step -> GetTrack() -> GetTrackStatus() != fAlive) { return; }

  G4LogicalVolume * preStepPointVolLog = 
       step -> GetPreStepPoint() -> GetPhysicalVolume() -> GetLogicalVolume();
  RegionInformation* preStepPointRegionInfo = (RegionInformation*) 
         (preStepPointVolLog -> GetRegion() -> GetUserInformation());

  G4LogicalVolume * postStepPointVolLog = 
       step -> GetPostStepPoint() -> GetPhysicalVolume() -> GetLogicalVolume();
  RegionInformation* postStepPointRegionInfo  = (RegionInformation*) 
         (postStepPointVolLog -> GetRegion() -> GetUserInformation());
 
  if(preStepPointRegionInfo -> IsWorld() && 
     postStepPointRegionInfo -> IsTarget()) { 
 
     G4Track* track = step -> GetTrack();

     G4String particle = track -> GetDefinition() -> GetParticleName();
     G4int trackID = track -> GetTrackID();
     G4double energy = track -> GetKineticEnergy(); 

     TrackInformation* trackInfo = 
                     (TrackInformation*) (track->GetUserInformation());
     
     if(trackInfo -> HitTarget()) {
        AnalysisManager::Instance() -> 
            ScoreExitingParticles(- trackInfo -> EscapeEnergy(), 0, 0, 0,
                                                           particle,trackID);  
     }
     else {
        AnalysisManager::Instance() -> ScoreEnteringParticles(energy, 0, 0, 0,
                                                              particle);  

     }
  }

  if(preStepPointRegionInfo -> IsTarget() && 
     postStepPointRegionInfo -> IsWorld()) { 

     G4Track* track = step -> GetTrack();

     G4String particle = track -> GetDefinition() -> GetParticleName();
     G4int trackID = track -> GetTrackID();
     G4double energy = track -> GetKineticEnergy(); 

     TrackInformation* trackInfo = 
                    (TrackInformation*) (track->GetUserInformation());
     
     trackInfo -> SetHitTarget();
     trackInfo -> SetEscapeEnergy(energy);

     AnalysisManager::Instance() -> ScoreExitingParticles(energy, 0, 0, 0,
                                                           particle,trackID);  
  }
}
