#ifndef G4ProtonAntiProtonReaction_h
#define G4ProtonAntiProtonReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4AntiProton.hh"

class G4ProtonAntiProtonReaction : public G4HadronicInteraction
{
  public: 
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theModel;
};

inline
G4VParticleChange * G4ProtonAntiProtonReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4AntiProton::AntiProton())
  {
    G4Exception("Calling G4ProtonAntiProtonReaction with particle other than p-bar!!!");
  }
  if(aTargetNucleus.GetZ() != 1)
  {
    G4Exception("Calling G4ProtonAntiProtonReaction for target other than Hydrogen!!!");
  }
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

#endif
