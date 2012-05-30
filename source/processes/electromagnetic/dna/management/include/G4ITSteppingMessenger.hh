#ifndef G4ITSTEPPINGMESSENGER_H
#define G4ITSTEPPINGMESSENGER_H

class G4ITStepManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ITSteppingMessenger: public G4UImessenger
{
  public:
    G4ITSteppingMessenger(G4ITStepManager* runMgr);
    ~G4ITSteppingMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    G4ITStepManager * fITStepManager;

  private: //commands
    G4UIdirectory*              fITDirectory;

    G4UIcmdWithADoubleAndUnit*  fEndTime;
    G4UIcmdWithAnInteger*       fVerboseCmd;
    G4UIcmdWithAnInteger*       fMaxStepNumber;
    G4UIcmdWithoutParameter*    fInitCmd;
    G4UIcmdWithoutParameter*    fProcessCmd;
    G4UIcmdWithAnInteger*       fMaxNULLTimeSteps;
};

#endif // G4ITSTEPPINGMESSENGER_H
