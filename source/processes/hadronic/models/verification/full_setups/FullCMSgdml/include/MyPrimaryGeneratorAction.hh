#ifndef MyPrimaryGeneratorAction_h
#define MyPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include <string>

class G4ParticleGun;
class G4Event;
class G4ParticleTable;


class MyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {

public:

  MyPrimaryGeneratorAction();
  ~MyPrimaryGeneratorAction();

public:

  void GeneratePrimaries( G4Event* anEvent );

private:

  G4ParticleGun* particleGun;
  G4ParticleTable* particleTable;

  static const int numberCandidateParticles;
  static const std::string nameParticlesVector[];

  const double infBeamEnergy;
  const double supBeamEnergy;

};

#endif


