#include "G4PScoreProcess.hh"
#include "G4VPScorer.hh"
#include "G4PStep.hh"
#include "G4VParallelStepper.hh"

G4PScoreProcess::G4PScoreProcess(G4VParallelStepper  &astepper,
				 G4VPScorer &aScorer,
				 const G4String &aName) :
  G4VProcess(aName), 
  fPstepper(astepper),
  fScorer(aScorer){
  G4VProcess::pParticleChange = new G4ParticleChange;
}

G4PScoreProcess::~G4PScoreProcess(){
  delete pParticleChange;
}

G4double G4PScoreProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition){
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
G4PScoreProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep){
  pParticleChange->Initialize(aTrack);
  
  fScorer.Score(aStep, fPstepper.GetPStep()); 
  
  return G4VProcess::pParticleChange;
}
