
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
