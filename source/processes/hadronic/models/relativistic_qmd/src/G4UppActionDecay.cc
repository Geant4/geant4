
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
  G4UppTrackChange* theChange = new G4UppTrackChange;
  theChange->oldParticles.push_back(particlePtr);
  G4KineticTrackVector* newPart = particlePtr->Decay();

  if (newPart==NULL) cout << "NULL!" << endl;
  // cout << "newpart: " << newPart->entries() << endl; // I hate rw
  for (G4int k=0; k<newPart->entries(); k++) {
    // cout << "DECPROD: " << (*newPart)[k]->GetDefinition()->GetParticleName() << endl;
    G4UppTrack* TrackPtr = new G4UppTrack(*(*newPart)[k]);
    TrackPtr->setLocalTime(allTracks.getGlobalTime());
    theChange->newParticles.push_back(TrackPtr);
  }

  return theChange;
}


G4bool G4UppActionDecay::isValid() const
{
  return !particlePtr->hasChanged();
}


void G4UppActionDecay::dump() const
{
  G4cout << "Action: DECAY ( ";
  G4cout << particlePtr->GetDefinition()->GetParticleName() << ",";
  G4cout << particlePtr->GetActualMass() << " ";
  G4cout << ") at " << getActionTime()*c_light/fermi << G4endl;
}
