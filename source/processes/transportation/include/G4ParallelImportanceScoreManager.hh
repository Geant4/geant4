#ifndef G4ParallelImportanceScoreManager_hh
#define G4ParallelImportanceScoreManager_hh

#include "G4ParallelImportanceManager.hh"
#include "G4ParallelWorld.hh"
#include "G4ProcessPlacer.hh"

class G4VPScorer;
class G4PScoreProcess;

class G4ParallelImportanceScoreManager :
  private G4ParallelImportanceManager {
public:
  G4ParallelImportanceScoreManager(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename);
  G4ParallelImportanceScoreManager(G4VIStore &iw, 
				   G4VPScorer &ascorer,
				   const G4String &particlename,
				   G4VImportanceAlgorithm &ialg);
  ~G4ParallelImportanceScoreManager();

  G4PScoreProcess *CreateParallelScoreProcess();
  void Initialize();

private:
  G4ParallelImportanceScoreManager(const G4ParallelImportanceScoreManager &);
  G4ParallelImportanceScoreManager &
  operator=(const G4ParallelImportanceScoreManager &);

  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScoreProcess;
};


#endif
