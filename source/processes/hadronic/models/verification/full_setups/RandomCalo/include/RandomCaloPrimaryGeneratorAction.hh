#ifndef RandomCaloPrimaryGeneratorAction_h
#define RandomCaloPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include <string>

class G4ParticleGun;
class G4Event;
class G4ParticleTable;


class RandomCaloPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  RandomCaloPrimaryGeneratorAction();
  ~RandomCaloPrimaryGeneratorAction();

public:

  void GeneratePrimaries( G4Event* anEvent );

private:

  G4ParticleGun* particleGun;
  G4ParticleTable* particleTable;

  static const int numberCandidateParticles;
  static const std::string nameParticlesVector[];

  const double infBeamEnergy;
  const double supBeamEnergy;

  const double rMax;
  const double zMax;
 
  const double energyThresholdNoBiasBelow;
  // For beam energies below this value no biasing is applied
  // to the event. If instead the beam energy is at least large
  // as this constant, then electrons/positrons/gammas/neutrons
  // are biased in the event. This biasing is applied by the 
  // class RandomCaloStackingAction. 

};

#endif


