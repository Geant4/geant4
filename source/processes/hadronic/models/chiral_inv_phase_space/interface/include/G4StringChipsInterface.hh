#ifndef G4StringChipsInterface_h
#define G4StringChipsInterface_h

#include "G4VIntraNuclearTransportModel.hh"

class G4StringChipsInterface : public G4VIntraNuclearTransportModel
{
  public:
    virtual G4VParticleChange* ApplyYourself(const G4Track& aTrack, 
                                             G4Nucleus& theNucleus);

    virtual G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
                                               G4V3DNucleus* theNucleus); 
};
#endif
