
#ifndef TstDrawVox01PrimaryGeneratorMessenger_h
#define TstDrawVox01PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TstDrawVox01PrimaryGeneratorAction;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

class TstDrawVox01PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    TstDrawVox01PrimaryGeneratorMessenger(TstDrawVox01PrimaryGeneratorAction * myPGA);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    TstDrawVox01PrimaryGeneratorAction * myPGAction;
    G4UIcmdWithoutParameter * standardGun;
    G4UIcmdWithAString * randomGun;
};

#endif

