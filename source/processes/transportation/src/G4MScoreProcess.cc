#include "G4MScoreProcess.hh"
#include "G4VPScorer.hh"
#include "G4PStep.hh"
#include "G4VParallelStepper.hh"

G4MScoreProcess::G4MScoreProcess(G4VPScorer &aScorer,
				 const G4String &aName) :
  G4VProcess(aName), 
  fScorer(aScorer){
  G4VProcess::pParticleChange = new G4ParticleChange;
}

G4MScoreProcess::~G4MScoreProcess(){
  delete pParticleChange;
}

G4double G4MScoreProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition){
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
G4MScoreProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep){
  pParticleChange->Initialize(aTrack);

  if (aStep.GetStepLength() > kCarTolerance) {
    G4StepPoint *prepoint = aStep.GetPreStepPoint();
    G4StepPoint *postpoint = aStep.GetPostStepPoint();
  

    G4PTouchableKey prekey(*(prepoint->GetPhysicalVolume()), 
			   prepoint->GetTouchable()->GetReplicaNumber());
    G4PTouchableKey postkey(*(postpoint->GetPhysicalVolume()), 
			    postpoint->GetTouchable()->GetReplicaNumber());

    G4PStep pstep(prekey, postkey);
    pstep.fCrossBoundary = false;
    
    if (prekey != postkey) pstep.fCrossBoundary = true;
  
    fScorer.Score(aStep, pstep); 
  }
  return G4VProcess::pParticleChange;
}
