#include "Tst33WeightChangeProcess.hh"
#include "Randomize.hh"

Tst33WeightChangeProcess::Tst33WeightChangeProcess()
 : 
  G4VProcess("Tst33WeightChangeProcess"),
  aParticleChange(new G4ParticleChange)
{
  G4VProcess::pParticleChange = aParticleChange;
}

Tst33WeightChangeProcess::~Tst33WeightChangeProcess()
{
  delete aParticleChange;
  aParticleChange = 0;
  G4VProcess::pParticleChange = 0;
}

G4double Tst33WeightChangeProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
Tst33WeightChangeProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &aStep)
{
  aParticleChange->Initialize(aTrack);
  G4double w = aTrack.GetWeight();

  G4double actionProb = 0.1;
  G4double lowestRelativeWeight = 0.1;
  G4double relativeWeightRange = 1 - lowestRelativeWeight;

  if (G4UniformRand() < actionProb) {
    G4double newWeight = w * (lowestRelativeWeight + 
			      relativeWeightRange * G4UniformRand());
    
    aParticleChange->SetWeightChange(newWeight);
  }
  return aParticleChange;
}

const G4String &Tst33WeightChangeProcess::GetName() const {
  return theProcessName;
}


G4double Tst33WeightChangeProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
				      G4double  ,
				      G4double  ,
				      G4double& ,
				      G4GPILSelection*) {
  return -1.0;
}

G4double Tst33WeightChangeProcess::
AtRestGetPhysicalInteractionLength(const G4Track&,
				   G4ForceCondition*) {
  return -1.0;
}

G4VParticleChange* Tst33WeightChangeProcess::AtRestDoIt(const G4Track&,
					       const G4Step&) {
  return 0;
}

G4VParticleChange* Tst33WeightChangeProcess::AlongStepDoIt(const G4Track&,
						  const G4Step&) {
  return 0;
}

