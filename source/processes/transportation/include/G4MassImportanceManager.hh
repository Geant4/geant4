#ifndef G4MassImportanceManager_hh
#define G4MassImportanceManager_hh G4MassImportanceManager_hh

#include "globals.hh"

class G4VIStore;
class G4MassImportanceProcess;
class G4VImportanceAlgorithm;
class G4VProcess;

class G4MassImportanceManager {
public:
  G4MassImportanceManager(G4VIStore &aIstore,
			  const G4String &particlename);

  G4MassImportanceManager(G4VIStore &aIstore,
			  const G4String &particlename,
			  const G4VImportanceAlgorithm &algorithm);
  
  ~G4MassImportanceManager();

  G4MassImportanceProcess *CreateMassImportanceProcess();
  void Initialize();

private:
  G4MassImportanceManager(const G4MassImportanceManager &);
  G4MassImportanceManager &operator=(const G4MassImportanceManager &);

  G4VIStore &fIStore;
  G4String fParticleName;
  const G4VImportanceAlgorithm &fAlgorithm;
  bool fCreatedAlgorithm;
  G4MassImportanceProcess *fMassImportanceProcess;
};

#endif
