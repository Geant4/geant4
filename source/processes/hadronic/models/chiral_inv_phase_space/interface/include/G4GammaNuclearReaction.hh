#ifndef G4GammaNuclearReaction_h
#define G4GammaNuclearReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4Gamma.hh"

class G4GammaNuclearReaction : public G4HadronicInteraction
{
  public: 
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theModel;
};

inline
G4VParticleChange * G4GammaNuclearReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4Gamma::GammaDefinition())
  {
    G4Exception("Called G4GammaNuclearReaction for particle other than gamma");
  }
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

#endif
