#ifndef StatAccepTestStackingAction_H
#define StatAccepTestStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include <fstream>

class G4Track;


class StatAccepTestStackingAction : public G4UserStackingAction {

  public:

  StatAccepTestStackingAction();
  virtual ~StatAccepTestStackingAction();

  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  // For gammas, electrons, and positrons it does the following
  // simple russian roulette biasing method:
  // if a particle of the above types is below a certain energy
  // threshold, and if its weight is below a certain threshold,
  // then once every twice the particle is killed, otherwise its
  // weight is doubled. 
  // In the case of neutron, a brutal approach is implemented:
  // any neutron below a certain energy threshold is killed.

  virtual void PrepareNewEvent();
  // Reset the counters for gammas, electrons, and positrons at the
  // beginning of each new events.

  private:

  G4int numberGammas1;
  G4int numberElectrons1;
  G4int numberPositrons1;
  G4int numberGammas2;
  G4int numberElectrons2;
  G4int numberPositrons2;
  G4int numberNeutrons;

  const G4double gammaEnergyThreshold1;
  const G4double electronEnergyThreshold1;
  const G4double positronEnergyThreshold1;
  // First (lower) energy threshold used for the russian 
  // roulette biasing: the biasing is applied only for 
  // energies below the first (higher) threshold,
  // with the first (higher) biasing factor. 
  // We assume:  energyThreshold1 <= energyThreshold2
  //             biasingFactor1   >= biasingFactor2
 
  const G4double gammaWeightThreshold1;
  const G4double electronWeightThreshold1;
  const G4double positronWeightThreshold1;
  // Weight threshold used for the russian roulette biasing:
  // the first biasing factor is applied only if the weight 
  // is below threshold1.

  const G4int gammaBiasingFactor1;
  const G4int electronBiasingFactor1;
  const G4int positronBiasingFactor1;
  // First (higher) biasing factor for the russian roulette 
  // biasing: one particle out of a number equal to the first 
  // biasing factor is kept, and weighted with such a factor, 
  // whereas the others are killed.

  const G4double gammaEnergyThreshold2;
  const G4double electronEnergyThreshold2;
  const G4double positronEnergyThreshold2;
  // Second (higher) energy threshold used for the russian 
  // roulette biasing: the biasing is applied only for energies 
  // below the second (higher) threshold, with the second 
  // (lower) biasing factor. 
  // We assume:  energyThreshold1 <= energyThreshold2
  //             biasingFactor1   >= biasingFactor2

  const G4double gammaWeightThreshold2;
  const G4double electronWeightThreshold2;
  const G4double positronWeightThreshold2;
  // Weight threshold used for the russian roulette biasing:
  // the second biasing factor is applied only if the weight 
  // is below the threshold2.

  const G4int gammaBiasingFactor2;
  const G4int electronBiasingFactor2;
  const G4int positronBiasingFactor2;
  // Second (lower) biasing factor for the russian roulette 
  // biasing: one particle out of a number equal to the second 
  // biasing factor is kept, and weighted with such a factor, 
  // whereas the others are killed.

  //ALB std::ofstream * filePtr;
  // Write out the weights of gammas, electrons, and positrons.

  enum { Nmax = 10 };

  G4int weightGammaVec[ Nmax ];
  G4int weightElectronVec[ Nmax ];
  G4int weightPositronVec[ Nmax ];
  // To see the distribution of weights: because the weights
  // can only be power of the first (lower) biasing factor 
  // (1, biasingFactor1^2, ...), 
  // it is enough to register the multiplicity for each exponent.

  const G4double neutronKillingEnergyThreshold;
  // Energy threshold below which a neutron is always killed
  // (so no biasing here, just a brutal cut!)

  const G4double neutronEnergyThreshold;
  const G4double neutronWeightThreshold;
  const G4int neutronBiasingFactor;
  // Parameters for a russian roulette biasing for neutrons.

  G4int weightNeutronVec[ Nmax ];
  // Distribution of weights for neutrons.

};


#endif

