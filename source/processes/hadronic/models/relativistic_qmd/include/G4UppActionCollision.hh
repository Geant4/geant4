
#ifndef G4UPPACTIONCOLLISION_H
#define G4UPPACTIONCOLLISION_H


#include "G4VUppAction.hh"
#include "G4VScatterer.hh"


class G4UppActionCollision : public G4VUppAction
{
public:

  G4UppActionCollision(const G4double collisionTime, 
		       const G4UppTrackVector& allTracks,
		       const G4UppTrackVector& collidingParticles,
		       G4VScatterer& aScatterer);

  G4bool isValid() const;

  G4UppTrackChange* perform(const G4UppTrackVector& allTracks) const;

  void dump() const;

private:

  const G4UppTrackVector* allTracksPtr;
  G4UppTrackVector collPart;
  G4VScatterer* scattererPtr;

};


#endif // G4UPPACTIONCOLLISION_H
