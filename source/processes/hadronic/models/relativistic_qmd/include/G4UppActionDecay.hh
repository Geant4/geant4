
#ifndef G4UPPACTIONDECAY_H
#define G4UPPACTIONDECAY_H


#include "G4VUppAction.hh"


class G4UppActionDecay : public G4VUppAction
{
public:

  G4UppActionDecay(const G4double decayTime, 
		   const G4UppTrack& decayingParticle);

  G4bool isValid() const;

  G4UppTrackChange* perform(const G4UppTrackVector& allTracks) const;

private:

  const G4UppTrack* particlePtr;

};


#endif // G4UPPACTIONDECAY_H
