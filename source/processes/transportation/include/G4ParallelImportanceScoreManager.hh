#ifndef G4ParallelImportanceScoreManager_hh
#define G4ParallelImportanceScoreManager_hh

#include "globals.hh"

class G4VIStore;
class G4VPScorer;
class G4VImportanceAlgorithm;
class G4PScoreProcess;
class G4ParallelManager;
class G4ParallelImportanceManager;

class G4ParallelImportanceScoreManager{
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

  G4ParallelManager &fParallelManager;
  G4ParallelImportanceManager &fParallelImportanceManager;
  G4VPScorer &fPScorer;
  G4PScoreProcess *fPScoreProcess;
};


#endif
