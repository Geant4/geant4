#ifndef G4MassScoreManager_hh
#define G4MassScoreManager_hh G4MassScoreManager_hh

#include "globals.hh"

class G4VProcess;
class G4VPScorer;
class G4MScoreProcess;

class G4MassScoreManager {
public:
  G4MassScoreManager(G4VPScorer &ascorer, const G4String &particlename);
  ~G4MassScoreManager();
  G4MScoreProcess *CreateMassScoreProcess();
  void Initialize();
private:
  G4MassScoreManager(const G4MassScoreManager &);
  G4MassScoreManager &operator=(const G4MassScoreManager &);
  G4VPScorer &fScorer;
  G4String fParticleName;
  G4MScoreProcess *fMScoreProcess;
};

#endif
