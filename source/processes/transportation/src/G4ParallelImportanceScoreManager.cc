#include "G4ParallelImportanceScoreManager.hh"
#include "G4VPScorer.hh"
#include "G4PScoreProcess.hh"

G4ParallelImportanceScoreManager::
G4ParallelImportanceScoreManager(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename):
  G4ParallelImportanceManager(is, particlename),
  fPScorer(ascorer),
  fPScoreProcess(0)
{}



G4ParallelImportanceScoreManager::
G4ParallelImportanceScoreManager(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename,
				 G4VImportanceAlgorithm &ialg):
  G4ParallelImportanceManager(is, particlename, ialg),
  fPScorer(ascorer),
  fPScoreProcess(0)
{}

G4ParallelImportanceScoreManager::
~G4ParallelImportanceScoreManager() {
  if (fPScoreProcess) delete  fPScoreProcess;
}

G4PScoreProcess *
G4ParallelImportanceScoreManager::CreateParallelScoreProcess(){
  if (!fPScoreProcess) {
    fPScoreProcess = 
      new G4PScoreProcess(G4ParallelManager::GetParallelWorld().
			  GetParallelStepper(), fPScorer);
  }
  return fPScoreProcess;
}

void G4ParallelImportanceScoreManager::Initialize() {
  G4ProcessPlacer placer(G4ParallelManager::GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  G4ParallelImportanceManager::Initialize();
}
