
#include "G4UppActionCollision.hh"


G4UppActionCollision::G4UppActionCollision(const G4double time,
					   const G4UppTrackVector& inParts)
{
  incomingPart = inParts;
  setActionTime(time);
}


G4int G4UppActionCollision::Perform(const G4UppTrackVector& t, 
				    G4UppInteraction& ia) const
{
  G4cout << "Collison of Particles " << G4endl;
  G4int arraySize = incomingPart.size();
  for (G4int i=0; i<arraySize; i++) {
    incomingPart[i]->clearLastPartner();
    incomingPart[i]->setChanged(true);
    for (G4int j=0; j<arraySize; j++) 
      incomingPart[i]->addLastPartner(incomingPart[j]);
  }
  cout << "collok" << endl;
  return 0;
}


G4bool G4UppActionCollision::isValid() const
{
  for (int i=0; i<incomingPart.size(); i++) 
    if (incomingPart[i]->hasChanged()) return false;
  return true;
}


void G4UppActionCollision::dump() const
{
  G4cout << "Action: COLL at " << getActionTime()/fermi << G4endl;
}


void G4UppActionCollision::dump(const G4UppTrackVector& v) const
{
  G4cout << "Action: COLL from ";
  G4cout << v.getIndex(incomingPart[0]) << " and ";
  G4cout << v.getIndex(incomingPart[1]);
  G4cout << " at " << getActionTime()/fermi << G4endl;
}

