#ifndef G4StringChipsInterface_h
#define G4StringChipsInterface_h

#include "G4VIntraNuclearTransportModel.hh"
#include "G4ChiralInvariantPhaseSpace.hh"

class G4StringChipsInterface : public G4VIntraNuclearTransportModel
{
  public:
    G4StringChipsInterface();
    virtual G4VParticleChange* ApplyYourself(const G4Track& aTrack, 
                                             G4Nucleus& theNucleus);

    virtual G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
                                               G4V3DNucleus* theNucleus); 
  private:
  
    G4ChiralInvariantPhaseSpace theModel;
    G4double theEnergyLossPerFermi;
};
#endif
