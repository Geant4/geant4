
#include "G4UppTrack.hh"


void G4UppTrack::dump() const
{
  G4cout << GetDefinition()->GetParticleName();
  G4cout << "  #Coll:" << NumberOfCollisions << G4endl;
}

G4bool G4UppTrack::isLastPartner(const G4UppTrack* PartPtr) const
{
  for (G4int i=0; i<lastPartners.size(); i++) {
    if (PartPtr==lastPartners[i]) return true;
  }
  return false;
}
