#ifndef G4DummyModel_h
#define G4DummyModel_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleChange.hh"

class G4DummyModel : public G4HadronicInteraction
{
  public:

  G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
  {
    return &theParticleChange;
  }

  private:

  G4ParticleChange theParticleChange;
  
};
#endif
