#ifndef G4VProcessPlacer_hh
#define G4VProcessPlacer_hh G4VProcessPlacer_hh
#include "globals.hh"

class G4VProcess;
 
class G4VProcessPlacer{
public:
  G4VProcessPlacer(const G4String &particlename){}
  virtual ~G4VProcessPlacer(){}
  virtual void AddProcessAsLastDoIt(G4VProcess *process) = 0;
  virtual void AddProcessAsSecondDoIt(G4VProcess *process) = 0;
};

#endif
