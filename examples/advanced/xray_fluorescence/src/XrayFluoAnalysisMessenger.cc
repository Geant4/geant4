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
  
  //   Commands for the 1D histograms 
     
  //   The Draw command gives the possibility to draw the 1d histograms 
  //   at every event.

  //   The Save command gives the possibility to save the 1d histograms in
  //   two separate PostScript files at the end of the run.
 

  Histo1DDrawCmd = new G4UIcmdWithAString("/analysis/histo1dDraw",this);
  Histo1DDrawCmd->SetGuidance("Enable the drawing of the 1d histograms every event.");
  Histo1DDrawCmd->SetGuidance("Choice: disable, enable(default)");
  Histo1DDrawCmd->SetParameterName("choice",true);
  Histo1DDrawCmd->SetDefaultValue("enable");
  Histo1DDrawCmd->SetCandidates("disable enable");
  Histo1DDrawCmd->AvailableForStates(Idle);

  Histo1DSaveCmd = new G4UIcmdWithAString("/analysis/histo1dSave",this);
  Histo1DSaveCmd->SetGuidance("Enable the saving of the 1d histograms every run.");
  Histo1DSaveCmd->SetGuidance("Choice: disable, enable(default)");
  Histo1DSaveCmd->SetParameterName("choice",true);
  Histo1DSaveCmd->SetDefaultValue("enable");
  Histo1DSaveCmd->SetCandidates("disable enable");
  Histo1DSaveCmd->AvailableForStates(Idle);

  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisMessenger::~XrayFluoAnalysisMessenger()
{
  delete Histo1DDrawCmd; 
  delete Histo1DSaveCmd; 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // 1D Histograms

  if( command == Histo1DDrawCmd )
    { XrayFluoAnalysis->SetHisto1DDraw(newValue);}

  if( command == Histo1DSaveCmd )
    { XrayFluoAnalysis->SetHisto1DSave(newValue);}
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif













