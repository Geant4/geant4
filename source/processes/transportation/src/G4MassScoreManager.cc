#include "G4MassScoreManager.hh"
#include "G4MScoreProcess.hh"
#include "G4ProcessPlacer.hh"


G4MassScoreManager::G4MassScoreManager(G4VPScorer &ascorer,
					const G4String &particlename):
  fScorer(ascorer),
  fParticleName(particlename),
  fMScoreProcess(0)
{}

G4MassScoreManager::~G4MassScoreManager(){
  if (fMScoreProcess) delete fMScoreProcess;
}



G4MScoreProcess *G4MassScoreManager::CreateMassScoreProcess(){
  if (!fMScoreProcess) {
    fMScoreProcess = new G4MScoreProcess(fScorer);
  }
  return fMScoreProcess;
}
void G4MassScoreManager::Initialize(){
  G4ProcessPlacer placer(fParticleName);
  placer.AddProcessAsSecondDoIt(CreateMassScoreProcess());
}





