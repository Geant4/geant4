// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4PImportanceWWindowScoreManager.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelManager.hh"
#include "G4PImportanceWWindowScoreManager.hh"
#include "G4ParallelImportanceManager.hh"
#include "G4ParallelWorld.hh"
#include "G4ProcessPlacer.hh"
#include "G4VIStore.hh"
#include "G4VPScorer.hh"
#include "G4PScoreProcess.hh"
#include "G4ParallelWeightWindowProcess.hh"
#include "G4WeightWindowAlgorithm.hh"

G4PImportanceWWindowScoreManager::
G4PImportanceWWindowScoreManager(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename,
				 G4VWeightWindowAlgorithm &wwalg)
 : fParallelManager(*(new G4ParallelManager(is.GetWorldVolume(), particlename))),
   fParallelImportanceManager(*(new G4ParallelImportanceManager(is, fParallelManager))),
   fPScorer(ascorer),
   fIstore(is),
   fPScoreProcess(0),
   fPWeightWindowProcess(0),
   fWWAlgorithm(wwalg)
{}

G4PImportanceWWindowScoreManager::
G4PImportanceWWindowScoreManager(G4VIStore &is, 
				 G4VPScorer &ascorer,
				 const G4String &particlename,
				 G4VWeightWindowAlgorithm &wwalg,
				 G4VImportanceAlgorithm &ialg)
 : fParallelManager(*(new G4ParallelManager(is.GetWorldVolume(), particlename))),
   fParallelImportanceManager(*(new 
			       G4ParallelImportanceManager(is,
							   ialg,
							   fParallelManager))),
   fPScorer(ascorer),
   fIstore(is),
   fPScoreProcess(0),
   fPWeightWindowProcess(0),
   fWWAlgorithm(wwalg)
{}

G4PImportanceWWindowScoreManager::~G4PImportanceWWindowScoreManager()
{
  if (fPScoreProcess) {
    G4ProcessPlacer placer(fParallelManager.GetParticleName());
    placer.RemoveProcess(fPScoreProcess);
    delete  fPScoreProcess;
  }
  if (fPWeightWindowProcess) {
    G4ProcessPlacer placer(fParallelManager.GetParticleName());
    placer.RemoveProcess(fPWeightWindowProcess);
    delete fPWeightWindowProcess;
  }
  delete &fParallelImportanceManager;
  delete &fParallelManager;
}

G4PScoreProcess *
G4PImportanceWWindowScoreManager::CreateParallelScoreProcess()
{
  if (!fPScoreProcess) {
    fPScoreProcess = 
      new G4PScoreProcess(fParallelManager.GetParallelWorld().
			  GetParallelStepper(), 
			  fPScorer);
  }
  return fPScoreProcess;
}

G4VProcess *
G4PImportanceWWindowScoreManager::CreateWeightWindowProcess()
{
  if (!fPWeightWindowProcess) {
    fPWeightWindowProcess = 
      new G4ParallelWeightWindowProcess(fIstore,
					fParallelManager.GetParallelWorld().
					GetParallelStepper(),
					fWWAlgorithm);
  }
  return fPWeightWindowProcess;
}

void G4PImportanceWWindowScoreManager::Initialize()
{
  G4ProcessPlacer placer(fParallelManager.GetParticleName());
  placer.AddProcessAsSecondDoIt(CreateParallelScoreProcess());
  fParallelImportanceManager.Initialize();
  placer.AddProcessAsSecondDoIt(CreateWeightWindowProcess());
}
