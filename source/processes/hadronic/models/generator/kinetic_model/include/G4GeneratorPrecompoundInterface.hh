
#ifndef G4GeneratorPrecompoundInterface_h
#define G4GeneratorPrecompoundInterface_h 1

#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"

class G4GeneratorPrecompoundInterface : public G4VIntraNuclearTransportModel 
{
public:
   G4GeneratorPrecompoundInterface(){}      
   ~G4GeneratorPrecompoundInterface(){}

private:
   G4int operator==(G4GeneratorPrecompoundInterface& right) {return (this == &right);}
   G4int operator!=(G4GeneratorPrecompoundInterface& right) {return (this != &right);}
      
public:
   G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus);
   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus);


private:   
};

#endif // G4GeneratorPrecompoundInterface_h


