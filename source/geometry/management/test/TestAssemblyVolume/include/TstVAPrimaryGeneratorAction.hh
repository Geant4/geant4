
#ifndef TstVAPrimaryGeneratorAction_h
#define TstVAPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4VisExtent.hh"

class G4ParticleGun;
class G4Event;
class TstVAPrimaryGeneratorMessenger;
class G4VPhysicalVolume;

class TstVAPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  enum Action {
    standardGun,
    randomDirectionGun,
    randomPositionGun,
    randomPositionAndDirectionGun
  };

  TstVAPrimaryGeneratorAction();
  ~TstVAPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  void SelectPrimaryGeneratorAction (Action);

private:
  Action generatorAction;
  G4ParticleGun* particleGun;
  TstVAPrimaryGeneratorMessenger* messenger;
  G4VPhysicalVolume* worldVolume;
  G4VisExtent worldExtent;
};

#endif
