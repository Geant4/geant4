#ifndef MLEventActionMessenger_h
#define MLEventActionMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

class MLEventAction;
////////////////////////////////////////////////////////////////////////////////
//
class MLEventActionMessenger: public G4UImessenger
{
  public:
    MLEventActionMessenger(MLEventAction*);
   ~MLEventActionMessenger();
    
    void SetNewValue (G4UIcommand*, G4String);
    
  private:
    MLEventAction         *eventAction;
    G4UIcmdWithADouble    *CpuCmd;
    G4UIcmdWithAString    *DrawCmd;
    G4UIcmdWithAnInteger  *PrintCmd;
    G4UIcmdWithAnInteger  *SeedCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif
