#include "G4ParallelImportanceManager.hh"
#include "G4ParallelManager.hh"
#include "G4ParallelWorld.hh"
#include "G4IStore.hh"
#include "G4ImportanceSampler.hh"
#include "G4ParallelImportanceProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ImportanceAlgorithm.hh"

G4ParallelImportanceManager::
G4ParallelImportanceManager(G4VIStore &is,
			    const G4String &particlename) :
  G4ParallelManager(is.GetWorldVolume(), particlename),
  fIalgorithm(*(new G4ImportanceAlgorithm())),
  fDeleteAlg(true),
  fSampler(new G4ImportanceSampler(fIalgorithm, 
				   G4ParallelManager::
				   GetParallelWorld().GetParallelStepper(),  
				   is)),
  fParallelImportanceProcess(0)
{}


  
G4ParallelImportanceManager::
G4ParallelImportanceManager(G4VIStore &is,
			    const G4String &particlename,
			    G4VImportanceAlgorithm &ialg): 
  G4ParallelManager(is.GetWorldVolume(), particlename),
  fIalgorithm(ialg),
  fDeleteAlg(false),
  fSampler(new G4ImportanceSampler(fIalgorithm, 
				   G4ParallelManager::
				   GetParallelWorld().GetParallelStepper(),  
				   is)),
  fParallelImportanceProcess(0)
{}
  
G4ParallelImportanceProcess *
G4ParallelImportanceManager::CreateParallelImportanceProcess(){
  if (!fParallelImportanceProcess) {
    fParallelImportanceProcess = 
      new G4ParallelImportanceProcess(*fSampler, 
				      G4ParallelManager::
				      GetParallelWorld().
				      GetGeoDriver(), 
				      G4ParallelManager::
				      GetParallelWorld().
				      GetParallelStepper());
  }
  return fParallelImportanceProcess;
}
void G4ParallelImportanceManager::Initialize(){
  G4ProcessPlacer placer(G4ParallelManager::GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelImportanceProcess());
}

G4ParallelImportanceManager::~G4ParallelImportanceManager(){
  if (fDeleteAlg) delete &fIalgorithm;
  if (fParallelImportanceProcess) delete fParallelImportanceProcess;
  delete fSampler;
}








