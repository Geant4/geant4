#ifndef G4CASCADEINTERFACE_H
#define G4CASCADEINTERFACE_H 1

//#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

// Class Description
// HETC/INUCL implementation of an intra-nuclear transport

class G4CascadeInterface : public G4VIntraNuclearTransportModel 
{
public:
   G4CascadeInterface(){}      
   ~G4CascadeInterface(){}

private:
   G4int operator==(G4CascadeInterface& right) {return (this == &right);}
   G4int operator!=(G4CascadeInterface& right) {return (this != &right);}
      
public:
   G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus);
   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);


private:   
};

#endif // G4CASCADEINTERFACE_H



