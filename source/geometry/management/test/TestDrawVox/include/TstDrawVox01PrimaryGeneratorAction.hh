
#ifndef TstDrawVox01PrimaryGeneratorAction_h
#define TstDrawVox01PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4VisExtent.hh"

class G4ParticleGun;
class G4Event;
class TstDrawVox01PrimaryGeneratorMessenger;
class G4VPhysicalVolume;

class TstDrawVox01PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:

  enum Action {
    standardGun,
    randomDirectionGun,
    randomPositionGun,
    randomPositionAndDirectionGun
  };

  TstDrawVox01PrimaryGeneratorAction();
  ~TstDrawVox01PrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  void SelectPrimaryGeneratorAction (Action);

private:
  Action generatorAction;
  G4ParticleGun* particleGun;
  TstDrawVox01PrimaryGeneratorMessenger* messenger;
  G4VPhysicalVolume* worldVolume;
  G4VisExtent worldExtent;
};

#endif
