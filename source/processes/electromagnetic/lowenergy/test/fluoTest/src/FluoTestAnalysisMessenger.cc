#ifdef G4ANALYSIS_USE
#include "FluoTestAnalysisMessenger.hh"

#include "FluoTestAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisMessenger::FluoTestAnalysisMessenger(FluoTestAnalysisManager* analysisManager)
  :FluoTestAnalysis(analysisManager)

{ 
  FluoTestAnalysisDir = new G4UIdirectory("/analysis/");
  FluoTestAnalysisDir->SetGuidance("esperimento analysis control.");
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisMessenger::~FluoTestAnalysisMessenger()
{
  
  delete FluoTestAnalysisDir; 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif













