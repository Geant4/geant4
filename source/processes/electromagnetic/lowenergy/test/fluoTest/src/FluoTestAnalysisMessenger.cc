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
  
  //   Commands for the 1D histograms (energy deposition in the last 
  //    layer and hits distribution along the sample)
     
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

  //    Commands for the 2D histograms (hits positions along the sample)
     
  //    The Draw command gives the possibility to draw the 1d histograms 
  //   at every event.

  //    The Save command gives the possibility to save the 1d histograms in
  //   two separate PostScript files at the end of the run.

  // Moreover there is the possibility to set the 2d histograms so
  //   that the info stored are true position ((x,z) or (y,z)
  //   coordinates with respect to the payload reference frame in mm) or
  //   the number of the Strip and the number of the Plane in which the
  //   hit occur. To note that this feature is just for visualization 
  //   purpouse since both the information are saved in the external ASCII
  //   file.
  
  //Histo2DDrawCmd = new G4UIcmdWithAString("/analysis/histo2dDraw",this);
  //Histo2DDrawCmd->SetGuidance("Enable the drawing of the 2d histograms every events.");
  // Histo2DDrawCmd->SetGuidance("Choice: disable, enable(default)");
  //Histo2DDrawCmd->SetParameterName("choice",true);
  //Histo2DDrawCmd->SetDefaultValue("enable");
  //Histo2DDrawCmd->SetCandidates("disable enable");
  //Histo2DDrawCmd->AvailableForStates(Idle);

  //Histo2DSaveCmd = new G4UIcmdWithAString("/analysis/histo2dSave",this);
  //Histo2DSaveCmd->SetGuidance("Enable the saving of the 2d histograms every run.");
  //Histo2DSaveCmd->SetGuidance("Choice: disable, enable(default)");
  //Histo2DSaveCmd->SetParameterName("choice",true);
  //Histo2DSaveCmd->SetDefaultValue("enable");
  //Histo2DSaveCmd->SetCandidates("disable enable");
  //Histo2DSaveCmd->AvailableForStates(Idle);

  //Histo2DModeCmd = new G4UIcmdWithAString("/analysis/histo2dMode",this);
  //Histo2DModeCmd->SetGuidance("Select the mode for the 2d histograms.");
  //Histo2DModeCmd->SetGuidance("Choice: position, strip(default)");
  //Histo2DModeCmd->SetGuidance("position -> the histo is filled with true positions in mm");
  //Histo2DModeCmd->SetGuidance("strip -> the histo is filled with the number of the strip and the plane");
  //Histo2DModeCmd->SetParameterName("choice",true);
  //Histo2DModeCmd->SetDefaultValue("strip");
  //Histo2DModeCmd->SetCandidates("position strip");
  //Histo2DModeCmd->AvailableForStates(Idle);
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestAnalysisMessenger::~FluoTestAnalysisMessenger()
{
  delete Histo1DDrawCmd; 
  delete Histo1DSaveCmd; 
  delete FluoTestAnalysisDir; //aggiunto
  //  delete Histo2DDrawCmd; 
  //delete Histo2DSaveCmd; 
  //delete Histo2DModeCmd; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // 1D Histograms

  if( command == Histo1DDrawCmd )
    { FluoTestAnalysis->SetHisto1DDraw(newValue);}

  if( command == Histo1DSaveCmd )
    { FluoTestAnalysis->SetHisto1DSave(newValue);}

  // 2D Histograms

  //if( command == Histo2DDrawCmd )
  // { FluoTestAnalysis->SetHisto2DDraw(newValue);}

  //if( command == Histo2DSaveCmd )
  //{ FluoTestAnalysis->SetHisto2DSave(newValue);}

  //if( command == Histo2DModeCmd )
  //  { FluoTestAnalysis->SetHisto2DMode(newValue);}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif













