#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisMessenger.hh"

#include "XrayFluoAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisMessenger::XrayFluoAnalysisMessenger(XrayFluoAnalysisManager* analysisManager)
  :XrayFluoAnalysis(analysisManager)

{ 
  XrayFluoAnalysisDir = new G4UIdirectory("/analysis/");
  XrayFluoAnalysisDir->SetGuidance("esperimento analysis control.");
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisMessenger::~XrayFluoAnalysisMessenger()
{
  
  delete XrayFluoAnalysisDir; 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif













