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

  inline G4int getNumberOfGammas();
  inline G4int getNumberOfElectrons();
  inline G4int getNumberOfPositrons();

  private:

  G4int numberGammas;
  G4int numberElectrons;
  G4int numberPositrons;

  const G4double gammaEnergyThreshold;
  const G4double electronEnergyThreshold;
  const G4double positronEnergyThreshold;
  // Energy threshold used for the russian roulette biasing:
  // the biasing is applied only for energies below the threshold.

  const G4double gammaWeightThreshold;
  const G4double electronWeightThreshold;
  const G4double positronWeightThreshold;
  // Weight threshold used for the russian roulette biasing:
  // the biasing is applied only if the weight is below the threshold.

  //ALB std::ofstream * filePtr;
  // Write out the weights of gammas, electrons, and positrons.

  enum { Nmax = 20 };
  G4int weightVec[ Nmax ];
  // To see the distribution of weights: because the weights
  // can only be power of 2 (1, 2, 4, 8, 16, ... 2^k), it is
  // enough to register the multiplicity for each exponent.

  const G4double neutronEnergyThreshold;
  // Energy threshold below which a neutron is killed.

};


inline G4int StatAccepTestStackingAction::getNumberOfGammas() {
  return numberGammas;
}

inline G4int StatAccepTestStackingAction::getNumberOfElectrons() {
  return numberElectrons;
}

inline G4int StatAccepTestStackingAction::getNumberOfPositrons() {
  return numberPositrons;
}


#endif

