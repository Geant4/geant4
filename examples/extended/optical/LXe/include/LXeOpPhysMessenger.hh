
#ifndef LXeOpPhysMessenger_h
#define LXeOpPhysMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class LXeOpticalPhysics;
class G4UIcmdWithADouble;

class LXeOpPhysMessenger: public G4UImessenger
{
public:
  LXeOpPhysMessenger(LXeOpticalPhysics*);
  ~LXeOpPhysMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  LXeOpticalPhysics* OpPhys;

  G4UIcmdWithADouble*          yieldCmd;

};

#endif

