
#include "G4UppActionCollision.hh"


G4UppActionCollision::G4UppActionCollision(const G4double collisionTime,
					   const G4UppTrackVector& allTracks,
					   const G4UppTrackVector& collidingParticles,
					   const G4VScatterer& aScatterer)
  : collPart(collidingParticles)
{
  setActionTime(collisionTime);
  allTracksPtr = &allTracks;
  // collPart = collidingParticles;
  scattererPtr = &aScatterer;
}


G4UppTrackChange* G4UppActionCollision::perform(const G4UppTrackVector& allTracks) const
{
  G4cout << "(debug) performing collison (";
  for (G4int j=0; j<collPart.size(); j++) {
    G4cout << collPart[j]->GetDefinition()->GetParticleName() << " ";
  }
  G4cout << ")" << G4endl;

  G4int arraySize = collPart.size();
  for (G4int i=0; i<arraySize; i++) {
    collPart[i]->clearLastInteractionPartners();
    collPart[i]->setChanged(true);
    for (G4int j=0; j<arraySize; j++) 
      collPart[i]->addLastInteractionPartner(collPart[j]);
  }
  // call scatterer here!
  return NULL;
}


G4bool G4UppActionCollision::isValid() const
{
  for (G4int i=0; i<collPart.size(); i++) 
    if (collPart[i]->hasChanged()) return false;
  return true;
}



void G4UppActionCollision::dump() const
{
  G4cout << "Action: COLL ( ";
  for (G4int i=0; i<collPart.size(); i++) {
    G4cout << allTracksPtr->getIndex(collPart[i]) << " ";
  }
  G4cout << ") at " << getActionTime()/fermi << G4endl;
}
