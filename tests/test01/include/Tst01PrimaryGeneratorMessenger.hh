
#ifndef Tst01PrimaryGeneratorMessenger_h
#define Tst01PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst01PrimaryGeneratorAction;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

class Tst01PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    Tst01PrimaryGeneratorMessenger(Tst01PrimaryGeneratorAction * myPGA);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst01PrimaryGeneratorAction * myPGAction;
    G4UIcmdWithoutParameter * standardGun;
    G4UIcmdWithAString * randomGun;
};

#endif

