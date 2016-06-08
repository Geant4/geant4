#ifndef G4ParametrizedHadronicVertex_h
#define G4ParametrizedHadronicVertex_h 1

#include "globals.hh"
#include "G4ParticleChange.hh"
#include "G4Nucleus.hh"
#include "G4ReactionProduct.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4Track.hh"

class G4ParametrizedHadronicVertex
{
  public:
   G4VParticleChange * ApplyYourself(G4Nucleus & theTarget, 
                                     const G4Track &thePhoton);

  private:
   G4LEPionPlusInelastic  theLowEPionPlus;
   G4LEPionMinusInelastic theLowEPionMinus;
   G4HEPionPlusInelastic  theHighEPionPlus;
   G4HEPionMinusInelastic theHighEPionMinus;
  
};
#endif
