 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"

extern ofstream outFile;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
XrayFluoRunAction::XrayFluoRunAction(XrayFluoAnalysisManager* aMgr)
  :analysisManager(aMgr)
{
}
#else
myRunAction::XrayFluoRunAction()
{
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoRunAction::~XrayFluoRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::BeginOfRunAction(const G4Run* aRun)
{
  char name[15];
  
  // Open the file for the tracks of this run
  sprintf(name,"Tracks_%d.dat", aRun->GetRunID());
  outFile.open(name);
  
  G4cout << "### Run " << aRun << " start." << G4endl;
  
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
#ifdef G4ANALYSIS_USE
  analysisManager->BeginOfRun();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoRunAction::EndOfRunAction(const G4Run* aRun )
{
  // Run ended, update the visualization
  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
  
  // Close the file with the hits information
  outFile.close();
  
  // If analysis is used, print out the histograms
#ifdef G4ANALYSIS_USE
  analysisManager->EndOfRun(aRun->GetRunID());
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







