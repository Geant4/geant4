
#include "G4UppActionDecay.hh"


G4UppActionDecay::G4UppActionDecay(const G4double decayTime,
				   const G4UppTrack& decayingParticle)
{
  particlePtr = &decayingParticle;
  setActionTime(decayTime);
}


G4UppTrackChange* G4UppActionDecay::perform(const G4UppTrackVector& allTracks) const
{
  G4cout << "(debug) performing decay of ";
  G4cout << particlePtr->GetDefinition()->GetParticleName() << G4endl;
  // insert code here
  return NULL;
}


G4bool G4UppActionDecay::isValid() const
{
  return !particlePtr->hasChanged();
}

