#ifndef G4BinaryLightIonReaction_h
#define G4BinaryLightIonReaction_h

#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4HadFinalState.hh"
#include "G4ExcitationHandler.hh"

class G4BinaryLightIonReaction : public G4HadronicInteraction 
{
  public:
    G4BinaryLightIonReaction();
    virtual ~G4BinaryLightIonReaction(){}
    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                              G4Nucleus& theNucleus);
  
  private:
    G4BinaryCascade theModel;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel theProjectileFragmentation;
    G4HadFinalState theResult;
    G4bool EnergyAndMomentumCorrector(G4ReactionProductVector* products,
    				G4LorentzVector& TotalCollisionMom);
};

#endif
