
#ifndef ExN04StackingActionMessenger_h
#define ExN04StackingActionMessenger_h 1

class ExN04StackingAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class ExN04StackingActionMessenger: public G4UImessenger
{
  public:
    ExN04StackingActionMessenger(ExN04StackingAction* msa);
    ~ExN04StackingActionMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    ExN04StackingAction * myAction;
    
  private: //commands
    G4UIcmdWithAnInteger * muonCmd;
    G4UIcmdWithAnInteger * isomuonCmd;
    G4UIcmdWithAnInteger * isoCmd;
    G4UIcmdWithADoubleAndUnit * roiCmd;
    
};

#endif


