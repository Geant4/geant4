
#ifndef Tst01PrimaryGeneratorAction_h
#define Tst01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4VisExtent.hh"

class G4ParticleGun;
class G4Event;
class Tst01PrimaryGeneratorMessenger;
class G4VPhysicalVolume;

class Tst01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  enum Action {
    standardGun,
    randomDirectionGun,
    randomPositionGun,
    randomPositionAndDirectionGun
  };

  Tst01PrimaryGeneratorAction();
  ~Tst01PrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  void SelectPrimaryGeneratorAction (Action);

private:
  Action generatorAction;
  G4ParticleGun* particleGun;
  Tst01PrimaryGeneratorMessenger* messenger;
  G4VPhysicalVolume* worldVolume;
  G4VisExtent worldExtent;
};

#endif
