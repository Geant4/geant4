

#include "G4UppInteraction.hh"
#include "G4KineticTrack.hh"


G4UppInteraction::G4UppInteraction(const G4double Time, const G4int i, const G4int j)
{
  IncomingParticleIndex.push_back(i);
  IncomingParticleIndex.push_back(j);
  InteractionTime = Time;
}


G4UppTrackChange G4UppInteraction::Perform()
{
  G4UppTrackChange aChange;
  return aChange;
}
