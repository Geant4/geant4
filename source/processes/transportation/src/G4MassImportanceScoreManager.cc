#include "G4MassImportanceScoreManager.hh"


G4MassImportanceScoreManager::
G4MassImportanceScoreManager(G4VIStore &aIstore,
			     G4VPScorer &ascorer,
			     const G4String &particlename):
  G4MassImportanceManager(aIstore, particlename),
  G4MassScoreManager(ascorer, particlename)
{}
  
G4MassImportanceScoreManager::
G4MassImportanceScoreManager(G4VIStore &aIstore,
			     G4VPScorer &ascorer,
			     const G4String &particlename,
			     const G4VImportanceAlgorithm &algorithm):
  G4MassImportanceManager(aIstore, particlename, algorithm),
  G4MassScoreManager(ascorer, particlename)
{}

void G4MassImportanceScoreManager::Initialize(){
  G4MassImportanceManager::Initialize();
  G4MassScoreManager::Initialize();
}
