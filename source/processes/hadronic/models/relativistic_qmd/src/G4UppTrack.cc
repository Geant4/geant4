
#include "G4UppTrack.hh"


void G4UppTrack::dump() const
{
  G4cout << GetDefinition()->GetParticleName();
  G4cout << "  nColl=" << numberOfCollisions << G4endl;
}


G4bool G4UppTrack::isLastInteractionPartner(const G4UppTrack* PartPtr) const
{
  for (G4int i=0; i<lastInteractionPartners.size(); i++) {
    if (PartPtr==lastInteractionPartners[i]) return true;
  }
  return false;
}
