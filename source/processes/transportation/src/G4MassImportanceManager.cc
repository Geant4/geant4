#include "G4MassImportanceManager.hh"
#include "G4MassImportanceProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ImportanceAlgorithm.hh"


G4MassImportanceManager::
G4MassImportanceManager(G4VIStore &aIstore,
			const G4String &particlename,
			const G4VImportanceAlgorithm &algorithm):
  fIStore(aIstore),
  fParticleName(particlename),
  fAlgorithm(algorithm),
  fCreatedAlgorithm(false),
  fMassImportanceProcess(0)
{}

G4MassImportanceManager::
G4MassImportanceManager(G4VIStore &aIstore,
			const G4String &particlename):
  fIStore(aIstore),
  fParticleName(particlename),
  fAlgorithm(*(new G4ImportanceAlgorithm)),
  fCreatedAlgorithm(true),
  fMassImportanceProcess(0)
{}


G4MassImportanceManager::~G4MassImportanceManager(){
  if (fMassImportanceProcess) delete fMassImportanceProcess;
  if (fCreatedAlgorithm) delete &fAlgorithm;
}


G4MassImportanceProcess *G4MassImportanceManager::
CreateMassImportanceProcess(){
  if (!fMassImportanceProcess) {
    fMassImportanceProcess =
      new G4MassImportanceProcess(fAlgorithm, fIStore);
  }
  return fMassImportanceProcess;
}

void G4MassImportanceManager::Initialize(){
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsLastDoIt(CreateMassImportanceProcess());
}
