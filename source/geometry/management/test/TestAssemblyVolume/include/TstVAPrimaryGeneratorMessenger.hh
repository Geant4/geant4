
#ifndef TstVAPrimaryGeneratorMessenger_h
#define TstVAPrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TstVAPrimaryGeneratorAction;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

class TstVAPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    TstVAPrimaryGeneratorMessenger(TstVAPrimaryGeneratorAction * myPGA);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    TstVAPrimaryGeneratorAction * myPGAction;
    G4UIcmdWithoutParameter * standardGun;
    G4UIcmdWithAString * randomGun;
};

#endif

