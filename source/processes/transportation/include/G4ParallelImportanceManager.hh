#ifndef G4ParallelImportanceManager_hh
#define G4ParallelImportanceManager_hh G4ParallelImportanceManager_hh

#include "globals.hh"
class G4ParallelManager;

class G4VIStore;
class G4VImportanceAlgorithm;
class G4VImportanceSampler;
class G4ParallelImportanceProcess;

class G4ParallelImportanceManager{
public:
  G4ParallelImportanceManager(G4VIStore &is, 
			      const G4String &particlename);
  G4ParallelImportanceManager(G4VIStore &is, 
			      G4ParallelManager &pmanager);
  G4ParallelImportanceManager(G4VIStore &is, 
			      const G4String &particlename,
			      G4VImportanceAlgorithm &ialg);
  G4ParallelImportanceManager(G4VIStore &is, 
			      G4VImportanceAlgorithm &ialg,
			      G4ParallelManager &pmanager);
  
  virtual ~G4ParallelImportanceManager();
  G4ParallelImportanceProcess *CreateParallelImportanceProcess();
  void Initialize();
private:
  G4ParallelImportanceManager(const G4ParallelImportanceManager &);
  G4ParallelImportanceManager &operator=(const G4ParallelImportanceManager &);

  G4ParallelManager &fParallelManager;
  bool fCreatedPM;
  G4VImportanceAlgorithm &fIalgorithm;
  G4bool fDeleteAlg;

  G4VImportanceSampler *fSampler;
  G4ParallelImportanceProcess *fParallelImportanceProcess;
};

#endif








