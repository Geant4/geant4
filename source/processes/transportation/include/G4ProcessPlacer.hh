#ifndef G4ProcessPlacer_hh
#define G4ProcessPlacer_hh
#include "globals.hh"
#include "G4VProcessPlacer.hh"

class G4ProcessManager;
class G4ProcessVector;

class G4ProcessPlacer : public G4VProcessPlacer{
public:
  G4ProcessPlacer(const G4String &particlename);
  void AddProcessAsLastDoIt(G4VProcess *process);
  void AddProcessAsSecondDoIt(G4VProcess *process);
private:
  G4ProcessManager &GetProcessManager();
  enum SecondOrLast {
    eSecond = 1,            
    eLast = 0
  };
  void AddProcessAs(G4VProcess *process, SecondOrLast);
  void PrintProcVec(G4ProcessVector* processVec);  

  G4String fParticleName;

};

#endif
