#include "G4ParallelManager.hh"
#include "G4ParallelImportanceScoreManager.hh"
#include "G4ParallelImportanceManager.hh"
#include "G4ParallelWorld.hh"
#include "G4ProcessPlacer.hh"
#include "G4VIStore.hh"
#include "G4VPScorer.hh"
#include "G4PScoreProcess.hh"

G4ParallelImportanceScoreManager::
G4ParallelImportanceScoreManager(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename):
  fParallelManager(*(new 
		     G4ParallelManager(is.GetWorldVolume(), particlename))),
  fParallelImportanceManager(*(new 
			       G4ParallelImportanceManager(is,
							   fParallelManager))),
  fPScorer(ascorer),
  fPScoreProcess(0)
{}



G4ParallelImportanceScoreManager::
G4ParallelImportanceScoreManager(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename,
				 G4VImportanceAlgorithm &ialg):
  fParallelManager(*(new 
		     G4ParallelManager(is.GetWorldVolume(), particlename))),
  fParallelImportanceManager(*(new 
			       G4ParallelImportanceManager(is,
							   ialg,
							   fParallelManager))),
  fPScorer(ascorer),
  fPScoreProcess(0)
{}

G4ParallelImportanceScoreManager::
~G4ParallelImportanceScoreManager() {
  delete &fParallelImportanceManager;
  delete &fParallelManager;
  if (fPScoreProcess) delete  fPScoreProcess;
}

G4PScoreProcess *
G4ParallelImportanceScoreManager::CreateParallelScoreProcess(){
  if (!fPScoreProcess) {
    fPScoreProcess = 
      new G4PScoreProcess(fParallelManager.GetParallelWorld().
			  GetParallelStepper(), fPScorer);
  }
  return fPScoreProcess;
}

void G4ParallelImportanceScoreManager::Initialize() {
  G4ProcessPlacer placer(fParallelManager.GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  fParallelImportanceManager.Initialize();
}
