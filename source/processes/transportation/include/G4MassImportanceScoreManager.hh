#ifndef G4MassImportanceScoreManager_hh
#define G4MassImportanceScoreManager_hh G4MassImportanceScoreManager_hh

#include "globals.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4MassImportanceManager;
class G4MassScoreManager;

class G4MassImportanceScoreManager {
public:
  G4MassImportanceScoreManager(G4VIStore &aIstore,
			       G4VPScorer &ascorer,
			       const G4String &particlename);
  
  G4MassImportanceScoreManager(G4VIStore &aIstore,
			       G4VPScorer &ascorer,
			       const G4String &particlename,
			       const G4VImportanceAlgorithm &algorithm);
  ~G4MassImportanceScoreManager();

  void Initialize();
private:
  G4MassImportanceScoreManager(const G4MassImportanceScoreManager &);
  G4MassImportanceScoreManager &
  operator=(const G4MassImportanceScoreManager &);

  G4MassImportanceManager *fMassImportanceManager;
  G4MassScoreManager *fMassScoreManager;
};

#endif
