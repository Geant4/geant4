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
  fParallelManager(*(new 
		     G4ParallelManager(is.GetWorldVolume(), particlename))),
  fCreatedPM(true),
  fIalgorithm(*(new G4ImportanceAlgorithm())),
  fDeleteAlg(true),
  fSampler(new G4ImportanceSampler(fIalgorithm, 
				   fParallelManager.
				   GetParallelWorld().GetParallelStepper(),  
				   is)),
  fParallelImportanceProcess(0)
{}

G4ParallelImportanceManager::
G4ParallelImportanceManager(G4VIStore &is, 
			    G4ParallelManager &pmanager):
  fParallelManager(pmanager),
  fCreatedPM(false),
  fIalgorithm(*(new G4ImportanceAlgorithm())),
  fDeleteAlg(true),
  fSampler(new G4ImportanceSampler(fIalgorithm, 
				   fParallelManager.
				   GetParallelWorld().GetParallelStepper(),  
				   is)),
  fParallelImportanceProcess(0)
{}


G4ParallelImportanceManager::
G4ParallelImportanceManager(G4VIStore &is,
			    const G4String &particlename,
			    G4VImportanceAlgorithm &ialg): 
  fParallelManager(*(new 
		     G4ParallelManager(is.GetWorldVolume(), particlename))),
  fCreatedPM(true),
  fIalgorithm(ialg),
  fDeleteAlg(false),
  fSampler(new G4ImportanceSampler(fIalgorithm, 
				   fParallelManager.
				   GetParallelWorld().GetParallelStepper(),  
				   is)),
  fParallelImportanceProcess(0)
{}

G4ParallelImportanceManager::
G4ParallelImportanceManager(G4VIStore &is, 
			    G4VImportanceAlgorithm &ialg,
			    G4ParallelManager &pmanager):
  fParallelManager(pmanager),
  fCreatedPM(false),
  fIalgorithm(ialg),
  fDeleteAlg(false),
  fSampler(new G4ImportanceSampler(fIalgorithm, 
				   fParallelManager.
				   GetParallelWorld().GetParallelStepper(),  
				   is)),
  fParallelImportanceProcess(0)
{}
  
  
  
G4ParallelImportanceProcess *
G4ParallelImportanceManager::CreateParallelImportanceProcess(){
  if (!fParallelImportanceProcess) {
    fParallelImportanceProcess = 
      new G4ParallelImportanceProcess(*fSampler, 
				      fParallelManager.
				      GetParallelWorld().
				      GetGeoDriver(), 
				      fParallelManager.
				      GetParallelWorld().
				      GetParallelStepper());
  }
  return fParallelImportanceProcess;
}
void G4ParallelImportanceManager::Initialize(){
  G4ProcessPlacer placer(fParallelManager.GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelImportanceProcess());
}

G4ParallelImportanceManager::~G4ParallelImportanceManager(){
  if (fCreatedPM) delete &fParallelManager;
  if (fDeleteAlg) delete &fIalgorithm;
  if (fParallelImportanceProcess) delete fParallelImportanceProcess;
  delete fSampler;
}








