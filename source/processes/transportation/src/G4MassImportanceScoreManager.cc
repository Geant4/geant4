#include "G4MassImportanceScoreManager.hh"
#include "G4MassImportanceManager.hh"
#include "G4MassScoreManager.hh"


G4MassImportanceScoreManager::
G4MassImportanceScoreManager(G4VIStore &aIstore,
			     G4VPScorer &ascorer,
			     const G4String &particlename):
  fMassImportanceManager(new G4MassImportanceManager(aIstore, particlename)),
  fMassScoreManager(new G4MassScoreManager(ascorer, particlename))
{}
  
G4MassImportanceScoreManager::
G4MassImportanceScoreManager(G4VIStore &aIstore,
			     G4VPScorer &ascorer,
			     const G4String &particlename,
			     const G4VImportanceAlgorithm &algorithm):
  fMassImportanceManager(new 
			 G4MassImportanceManager(aIstore, 
						 particlename, algorithm)),
  fMassScoreManager(new G4MassScoreManager(ascorer, particlename))
{}

G4MassImportanceScoreManager::~G4MassImportanceScoreManager(){
  delete fMassScoreManager;
  delete fMassImportanceManager;
}

void G4MassImportanceScoreManager::Initialize(){
  fMassScoreManager->Initialize();
  fMassImportanceManager->Initialize();
}
