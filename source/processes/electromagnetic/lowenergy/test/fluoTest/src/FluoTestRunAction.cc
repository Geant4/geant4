 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestRunAction.hh"
#include <stdlib.h>
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "FluoTestDataSet.hh"
#include "FluoTestNormalization.hh"
//#include "FluoTestResponse.hh"
#include "G4LogLogInterpolation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
FluoTestRunAction::FluoTestRunAction(FluoTestAnalysisManager* aMgr)
  :analysisManager(aMgr)
{
 
   
}
FluoTestRunAction::FluoTestRunAction()
{
      G4double min = 0.*MeV;
      G4double max = 0.10000E+06 *MeV;
      G4int nBins = 835;

      FluoTestNormalization* normalization = new FluoTestNormalization();
     
      dataSet = normalization->Normalize(min, max, nBins,"mercury_flx_solmin");
      //dataSet = normalization->Normalize(min, max, nBins,"B_flare");
      
      G4String fileName = "efficienza";
      G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();
      efficiencySet = new FluoTestDataSet(1,fileName,interpolation,keV,1);

 delete normalization;

}
#else
FluoTestRunAction::FluoTestRunAction()
{
}
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestRunAction::~FluoTestRunAction()
{
  
  //delete responseFunction;
  //responseFunction = 0;
  // delete dataSet;
  //dataSet = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestRunAction::BeginOfRunAction(const G4Run* aRun)
{
  
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

void FluoTestRunAction::EndOfRunAction(const G4Run* aRun )
{
  // Run ended, update the visualization
if (G4VVisManager::GetConcreteInstance()) {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }
 
 
// If analysis is used, print out the histograms
#ifdef G4ANALYSIS_USE
  analysisManager->EndOfRun(aRun->GetRunID());
#endif
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


const FluoTestDataSet* FluoTestRunAction::GetSet()
{
  return  dataSet;
}
const FluoTestDataSet* FluoTestRunAction::GetEfficiencySet()
{
  return efficiencySet;
}




