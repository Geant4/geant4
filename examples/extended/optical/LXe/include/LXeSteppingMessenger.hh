
#ifndef LXeSteppingMessenger_h
#define LXeSteppingMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class LXeSteppingAction;
class G4UIcmdWithABool;

class LXeSteppingMessenger: public G4UImessenger
{
public:
  LXeSteppingMessenger(LXeSteppingAction*);
  ~LXeSteppingMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  LXeSteppingAction*        stepping;
  G4UIcmdWithABool*  oneStepPrimariesCmd;
  
};

#endif

