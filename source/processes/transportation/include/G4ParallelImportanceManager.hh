#ifndef G4ParallelImportanceManager_hh
#define G4ParallelImportanceManager_hh G4ParallelImportanceManager_hh

#include "globals.hh"
#include "G4ParallelManager.hh"

class G4VIStore;
class G4VImportanceAlgorithm;
class G4VImportanceSampler;
class G4ParallelImportanceProcess;

class G4ParallelImportanceManager : private G4ParallelManager{
public:
  G4ParallelImportanceManager(G4VIStore &iw, 
			      const G4String &particlename);
  G4ParallelImportanceManager(G4VIStore &iw, 
			      const G4String &particlename,
			      G4VImportanceAlgorithm &ialg);
  virtual ~G4ParallelImportanceManager();
  G4ParallelImportanceProcess *CreateParallelImportanceProcess();
  void Initialize();
private:
  G4ParallelImportanceManager(const G4ParallelImportanceManager &);
  G4ParallelImportanceManager &operator=(const G4ParallelImportanceManager &);

  G4VImportanceAlgorithm &fIalgorithm;
  G4bool fDeleteAlg;

  G4VImportanceSampler *fSampler;
  G4ParallelImportanceProcess *fParallelImportanceProcess;
};

#endif

