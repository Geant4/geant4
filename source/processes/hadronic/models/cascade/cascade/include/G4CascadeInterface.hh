#ifndef G4CASCADEINTERFACE_H
#define G4CASCADEINTERFACE_H 1

// CLASS DESCRIPTION
// G4CascadeInterface defines an interface to HETC and INUCL 
// models of an medium energy (~ 0.5 - 5 GeV) intra-nuclear transport.
// If you have any questions, please contact 
// package writer aatos.heikkinen@cern.ch.

#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

class G4CascadeInterface : public G4VIntraNuclearTransportModel {

public:

  G4CascadeInterface();

  ~G4CascadeInterface(){
  }

  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);

  G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus); // Don't use this

private:

  G4int operator==(G4CascadeInterface& right) {
    return (this == &right);
  }

  G4int operator!=(G4CascadeInterface& right) {
    return (this != &right);
  }

  G4int verboseLevel;

};

#endif // G4CASCADEINTERFACE_H
