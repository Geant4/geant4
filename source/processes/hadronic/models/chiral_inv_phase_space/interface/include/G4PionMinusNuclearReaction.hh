#ifndef G4PionMinusNuclearReaction_h
#define G4PionMinusNuclearReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4PionMinus.hh"

class G4PionMinusNuclearReaction : public G4HadronicInteraction
{
  public: 
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theModel;
};

inline
G4VParticleChange * G4PionMinusNuclearReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4PionMinus::PionMinusDefinition())
  {
    G4Exception("Called G4PionMinusNuclearReaction for particle other than PionMinus");
  }
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

#endif
