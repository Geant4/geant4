#ifndef G4MassImportanceScoreManager_hh
#define G4MassImportanceScoreManager_hh G4MassImportanceScoreManager_hh

#include "G4MassImportanceManager.hh"
#include "G4MassScoreManager.hh"

class G4MassImportanceScoreManager : 
  private G4MassImportanceManager,
  private G4MassScoreManager {
public:
  G4MassImportanceScoreManager(G4VIStore &aIstore,
			       G4VPScorer &ascorer,
			       const G4String &particlename);
  
  G4MassImportanceScoreManager(G4VIStore &aIstore,
			       G4VPScorer &ascorer,
			       const G4String &particlename,
			       const G4VImportanceAlgorithm &algorithm);
  
  void Initialize();
private:
  G4MassImportanceScoreManager(const G4MassImportanceScoreManager &);
  G4MassImportanceScoreManager &
  operator=(const G4MassImportanceScoreManager &);

};

#endif
