
#ifndef LXeEventMessenger_h
#define LXeEventMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class LXeEventAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

class LXeEventMessenger: public G4UImessenger
{
public:
  LXeEventMessenger(LXeEventAction*);
  ~LXeEventMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  LXeEventAction*        LXeEvent;
  G4UIcmdWithAnInteger*  saveThresholdCmd;
  G4UIcmdWithAnInteger*  verboseCmd;
  G4UIcmdWithAnInteger*  pmtThresholdCmd;
  G4UIcmdWithABool*      forceDrawPhotonsCmd;
  G4UIcmdWithABool*      forceDrawNoPhotonsCmd;
};

#endif

