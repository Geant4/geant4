#ifndef G4BinaryLightIonReaction_h
#define G4BinaryLightIonReaction_h

#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4ParticleChange.hh"
#include "G4ExcitationHandler.hh"

class G4BinaryLightIonReaction : public G4HadronicInteraction 
{
  public:
    G4BinaryLightIonReaction();
    virtual ~G4BinaryLightIonReaction(){}
    virtual G4VParticleChange *ApplyYourself(const G4Track &aTrack, G4Nucleus & targetNucleus );
  
  private:
    G4BinaryCascade theModel;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel theProjectileFragmentation;
    G4ParticleChange theResult;
};

#endif
