#ifndef G4BinaryLightIonReaction_h
#define G4BinaryLightIonReaction_h

#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"

class G4BinaryLightIonReaction : public G4HadronicInteraction 
{
  public:
    virtual G4VParticleChange *ApplyYourself(const G4Track &aTrack, G4Nucleus & targetNucleus );
  
  private:
    G4BinaryCascade theModel;
    G4PreCompoundModel theProjectileFragmentation;
};

#endif
