#ifndef G4StringChipsParticleLevelInterface_h
#define G4StringChipsParticleLevelInterface_h

#include "G4VIntraNuclearTransportModel.hh"
#include "G4ChiralInvariantPhaseSpace.hh"

class G4StringChipsParticleLevelInterface : public G4VIntraNuclearTransportModel
{
  public:
    G4StringChipsParticleLevelInterface();
    virtual G4VParticleChange* ApplyYourself(const G4Track& aTrack, 
                                             G4Nucleus& theNucleus);

    virtual G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
                                               G4V3DNucleus* theNucleus); 
  private:
  
    G4ChiralInvariantPhaseSpace theModel;
    G4double theEnergyLossPerFermi;
    
    G4double fractionOfSingleQuasiFreeNucleons;
    G4double fractionOfPairedQuasiFreeNucleons;
    G4double clusteringCoefficient;
    G4double temperature;
    G4double halfTheStrangenessOfSee;
    G4double etaToEtaPrime;
    G4double fusionToExchange;
    G4int nop;
};
#endif
