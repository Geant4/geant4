
#ifndef ExN04PrimaryGeneratorAction_h
#define ExN04PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4VPrimaryGenerator;
class G4Event;
class ExN04PrimaryGeneratorMessenger;

class ExN04PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    ExN04PrimaryGeneratorAction();
    ~ExN04PrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event* anEvent);

  private:
    G4VPrimaryGenerator* HEPEvt;
    G4VPrimaryGenerator* particleGun;
    ExN04PrimaryGeneratorMessenger* messenger;
    G4bool useHEPEvt;

  public:
    inline void SetHEPEvtGenerator(G4bool f)
    { useHEPEvt = f; }
    inline G4bool GetHEPEvtGenerator()
    { return useHEPEvt; }
};

#endif


