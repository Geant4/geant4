#include "G4ParallelScoreManager.hh"
#include "G4ParallelManager.hh"
#include "G4PScoreProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ParallelWorld.hh"

G4ParallelScoreManager::
G4ParallelScoreManager(G4VPhysicalVolume &worldvolume,
		       const G4String &particlename,
		       G4VPScorer &scorer) :
  fParallelManager(*(new G4ParallelManager(worldvolume, particlename))),
  fPScorer(scorer),
  fPScorerProcess(0)
{}

G4ParallelScoreManager::
~G4ParallelScoreManager(){
  delete &fParallelManager;
  if (fPScorerProcess) delete  fPScorerProcess;
}


G4PScoreProcess *G4ParallelScoreManager::CreateParallelScoreProcess(){
  if (!fPScorerProcess) {
    fPScorerProcess = 
      new G4PScoreProcess(fParallelManager.GetParallelWorld().
			  GetParallelStepper(), fPScorer);
  }
  return fPScorerProcess;
}

void G4ParallelScoreManager::Initialize(){
  G4ProcessPlacer placer(fParallelManager.GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  fParallelManager.Initialize();
}

