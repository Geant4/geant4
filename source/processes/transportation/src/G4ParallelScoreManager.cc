#include "G4ParallelScoreManager.hh"
#include "G4ParallelManager.hh"
#include "G4PScoreProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ParallelWorld.hh"

G4ParallelScoreManager::
G4ParallelScoreManager(G4VPhysicalVolume &worldvolume,
		       const G4String &particlename,
		       G4VPScorer &scorer) :
  G4ParallelManager(worldvolume, particlename),
  fPScorer(scorer),
  fPScorerProcess(0)
{}

G4ParallelScoreManager::
~G4ParallelScoreManager(){
  if (fPScorerProcess) delete  fPScorerProcess;
}


G4PScoreProcess *G4ParallelScoreManager::CreateParallelScoreProcess(){
  if (!fPScorerProcess) {
    fPScorerProcess = 
      new G4PScoreProcess(G4ParallelManager::GetParallelWorld().
			  GetParallelStepper(), fPScorer);
  }
  return fPScorerProcess;
}

void G4ParallelScoreManager::Initialize(){
  G4ProcessPlacer placer(G4ParallelManager::GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  G4ParallelManager::Initialize();
}

