
#ifndef G4DymmyINC_h
#define G4DymmyINC_h 1

#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticleVector.hh"

class G4DymmyINC : public G4VIntraNuclearTransportModel 
{
public:
   G4DymmyINC(){}      
   G4DymmyINC(G4double anEnergy)
   {
     theEnergy = anEnergy;
   }      
   ~G4DymmyINC(){}

private:
   G4int operator==(G4DymmyINC& right) {return (this == &right);}
   G4int operator!=(G4DymmyINC& right) {return (this != &right);}
   
   G4double theEnergy;
      
public:
   G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus)
   { return new G4ParticleChange;}

   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
   {
     G4ReactionProductVector * theFinalResult = new G4ReactionProductVector(0);
     theSecondaries->clearAndDestroy();
     delete theSecondaries;
     return theFinalResult;
   }


private:   
};

#endif // G4DymmyINC_h


