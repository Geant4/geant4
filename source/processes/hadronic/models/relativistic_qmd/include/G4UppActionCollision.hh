
#ifndef G4UPPACTIONCOLLISION_H
#define G4UPPACTIONCOLLISION_H


#include "G4VUppAction.hh"


class G4UppActionCollision : public G4VUppAction
{
public:

  G4UppActionCollision(const G4double time, 
		       const G4UppTrackVector& inPart);
  G4bool isValid() const;
  G4int Perform(const G4UppTrackVector& t) const { return 1; }
  G4int Perform(const G4UppTrackVector& t, G4UppInteraction& i) const;
  void dump() const;
  void dump(const G4UppTrackVector& t) const;

private:

  G4UppTrackVector incomingPart;

};


#endif // G4UPPACTIONCOLLISION_H
