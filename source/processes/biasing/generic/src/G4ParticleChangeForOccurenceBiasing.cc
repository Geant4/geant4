#include "G4ParticleChangeForOccurenceBiasing.hh"

G4ParticleChangeForOccurenceBiasing::G4ParticleChangeForOccurenceBiasing(G4String name)
  : G4VParticleChange(),
    fName(name),
    fWrappedParticleChange(0),
    fOccurenceWeightForNonInteraction(-1.0),
    fOccurenceWeightForInteraction(-1.0)
{}

G4ParticleChangeForOccurenceBiasing::~G4ParticleChangeForOccurenceBiasing()
{}

void G4ParticleChangeForOccurenceBiasing::SetWrappedParticleChange(G4VParticleChange* wpc)
{
  fWrappedParticleChange = wpc;
}

void G4ParticleChangeForOccurenceBiasing::StealSecondaries()
{
  SetNumberOfSecondaries( fWrappedParticleChange->GetNumberOfSecondaries() );
  for (G4int isecond = 0; isecond < fWrappedParticleChange->GetNumberOfSecondaries(); isecond++) 
    {
      G4Track* secondary =  fWrappedParticleChange->GetSecondary(isecond);
      secondary->SetWeight ( secondary->GetWeight() * fOccurenceWeightForInteraction );
      AddSecondary( secondary );
    }
  fWrappedParticleChange->Clear();
}


G4Step* G4ParticleChangeForOccurenceBiasing::UpdateStepForAtRest(G4Step* step)
{
  return step;
}

G4Step* G4ParticleChangeForOccurenceBiasing::UpdateStepForAlongStep(G4Step* step)
{
  // -- make particle change of wrapped process to apply its changes:
  if ( fWrappedParticleChange ) fWrappedParticleChange->UpdateStepForAlongStep( step );
  
  // -- multiply parent weight by weight due to occurence biasing:
  G4StepPoint* postStepPoint = step->GetPostStepPoint();
  //  G4cout << " part change along : w = " <<  postStepPoint->GetWeight() << " , wOcc =  " << fOccurenceWeightForNonInteraction << G4endl; // !
  postStepPoint->SetWeight( postStepPoint->GetWeight() * fOccurenceWeightForNonInteraction );
  
  return step;
}

G4Step* G4ParticleChangeForOccurenceBiasing::UpdateStepForPostStep(G4Step* step)
{
  // -- let make first wrapped process to apply its changes:
  fWrappedParticleChange->UpdateStepForPostStep(step);
  // -- then apply weight correction due to occurence biasing:
  G4StepPoint* postStepPoint = step->GetPostStepPoint();
  postStepPoint->SetWeight( postStepPoint->GetWeight() * fOccurenceWeightForInteraction );

  return step;
}
