
#ifndef ExN04PrimaryGeneratorMessenger_h
#define ExN04PrimaryGeneratorMessenger_h 1

class ExN04PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;

#include "G4UImessenger.hh"
#include "globals.hh"

class ExN04PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    ExN04PrimaryGeneratorMessenger(ExN04PrimaryGeneratorAction* mpga);
    ~ExN04PrimaryGeneratorMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    ExN04PrimaryGeneratorAction * myAction;
    
  private: //commands
    G4UIdirectory *             mydetDirectory;
    G4UIcmdWithAString *        genCmd;
    
};

#endif


