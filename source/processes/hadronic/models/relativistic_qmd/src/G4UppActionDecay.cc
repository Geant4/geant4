
#include "G4UppActionDecay.hh"


G4UppActionDecay::G4UppActionDecay(const G4double time,
				   G4UppTrack* pPtr)
{
  particlePtr = pPtr;
  setActionTime(time);
}


G4int G4UppActionDecay::Perform(const G4UppTrackVector& t, G4UppInteraction& i) const
{
  G4cout << "Decay of Particle " << particlePtr->GetDefinition()->GetParticleName();
  G4cout << G4endl;
  return 0;
}


G4bool G4UppActionDecay::isValid() const
{
  return !particlePtr->hasChanged();
}

