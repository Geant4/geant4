#ifndef PARN02JOB_HH
#define PARN02JOB_HH
#include "G4VtbbJob.hh"

class G4tbbRunManager;
class ExN02DetectorConstruction;

class ParN02Job : public G4VtbbJob {
public:
  ParN02Job(const G4String& macroFile);
  ~ParN02Job();
  //Interface implementation, from base class
  void CreateDetector(G4tbbRunManager* rm);
  void UserActions(G4tbbRunManager* rm);
  void InitSetup(G4tbbRunManager* rm );
  void JobPrepare(G4tbbRunManager* rm );
private:
  ExN02DetectorConstruction* detector;
};

#endif
