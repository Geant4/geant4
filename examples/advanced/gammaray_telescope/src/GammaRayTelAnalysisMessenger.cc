// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: GammaRayTelAnalysisMessenger.cc,v 1.1 2000-12-06 16:53:13 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      ------------ GammaRayTelAnalysisMessenger ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dic 2000)
//
// ************************************************************
#ifdef G4ANALYSIS_USE
#include "GammaRayTelAnalysisMessenger.hh"

#include "GammaRayTelAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysisMessenger::GammaRayTelAnalysisMessenger(GammaRayTelAnalysisManager* analysisManager)
  :GammaRayTelAnalysis(analysisManager)

{ 
  GammaRayTelAnalysisDir = new G4UIdirectory("/analysis/");
  GammaRayTelAnalysisDir->SetGuidance("GammaRayTel analysis control.");
  
  /* 
     Commands for the 1D histograms (energy deposition in the last 
     TKR layer and hits distribution along the TKR)
     
     The Draw command gives the possibility to draw the 1d histograms 
     at every event.

     The Save command gives the possibility to save the 1d histograms in
     two separate PostScript files at the end of the run.
  */

  Histo1DDrawCmd = new G4UIcmdWithAString("/analysis/histo1dDraw",this);
  Histo1DDrawCmd->SetGuidance("Enable the drawing of the 1d histograms every event.");
  Histo1DDrawCmd->SetGuidance("Choice: disable, enable(default)");
  Histo1DDrawCmd->SetParameterName("choice",true);
  Histo1DDrawCmd->SetDefaultValue("ebable");
  Histo1DDrawCmd->SetCandidates("disable enable");
  Histo1DDrawCmd->AvailableForStates(Idle);

  Histo1DSaveCmd = new G4UIcmdWithAString("/analysis/histo1dSave",this);
  Histo1DSaveCmd->SetGuidance("Enable the saving of the 1d histograms every run.");
  Histo1DSaveCmd->SetGuidance("Choice: disable, enable(default)");
  Histo1DSaveCmd->SetParameterName("choice",true);
  Histo1DSaveCmd->SetDefaultValue("enable");
  Histo1DSaveCmd->SetCandidates("disable enable");
  Histo1DSaveCmd->AvailableForStates(Idle);

  /* 
     Commands for the 2D histograms (hits positions along the TKR)
     
     The Draw command gives the possibility to draw the 1d histograms 
     at every event.

     The Save command gives the possibility to save the 1d histograms in
     two separate PostScript files at the end of the run.
     
     Moreover there is the possibility to set the 2d histograms so
     that the info stored are true position ((x,z) or (y,z)
     coordinates with respect to the payload reference frame in mm) or
     the number of the Strip and the number of the Plane in which the
     hit occur. To note that this feature is just for visualization 
     purpouse since both the information are saved in the external ASCII
     file.
  */
  Histo2DDrawCmd = new G4UIcmdWithAString("/analysis/histo2dDraw",this);
  Histo2DDrawCmd->SetGuidance("Enable the drawing of the 2d histograms every events.");
  Histo2DDrawCmd->SetGuidance("Choice: disable, enable(default)");
  Histo2DDrawCmd->SetParameterName("choice",true);
  Histo2DDrawCmd->SetDefaultValue("enable");
  Histo2DDrawCmd->SetCandidates("disable enable");
  Histo2DDrawCmd->AvailableForStates(Idle);

  Histo2DSaveCmd = new G4UIcmdWithAString("/analysis/histo2dSave",this);
  Histo2DSaveCmd->SetGuidance("Enable the saving of the 2d histograms every run.");
  Histo2DSaveCmd->SetGuidance("Choice: disable, enable(default)");
  Histo2DSaveCmd->SetParameterName("choice",true);
  Histo2DSaveCmd->SetDefaultValue("enable");
  Histo2DSaveCmd->SetCandidates("disable enable");
  Histo2DSaveCmd->AvailableForStates(Idle);

  Histo2DModeCmd = new G4UIcmdWithAString("/analysis/histo2dMode",this);
  Histo2DModeCmd->SetGuidance("Select the mode for the 2d histograms.");
  Histo2DModeCmd->SetGuidance("Choice: position, strip(default)");
  Histo2DModeCmd->SetGuidance("position -> the histo is filled with true positions in mm");
  Histo2DModeCmd->SetGuidance("strip -> the histo is filled with the number of the strip and the plane");
  Histo2DModeCmd->SetParameterName("choice",true);
  Histo2DModeCmd->SetDefaultValue("strip");
  Histo2DModeCmd->SetCandidates("position strip");
  Histo2DModeCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysisMessenger::~GammaRayTelAnalysisMessenger()
{
  delete Histo1DDrawCmd; 
  delete Histo1DSaveCmd; 
  delete Histo2DDrawCmd; 
  delete Histo2DSaveCmd; 
  delete Histo2DModeCmd; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  // 1D Histograms

  if( command == Histo1DDrawCmd )
    { GammaRayTelAnalysis->SetHisto1DDraw(newValue);}

  if( command == Histo1DSaveCmd )
    { GammaRayTelAnalysis->SetHisto1DSave(newValue);}

  // 2D Histograms

  if( command == Histo2DDrawCmd )
    { GammaRayTelAnalysis->SetHisto2DDraw(newValue);}

  if( command == Histo2DSaveCmd )
    { GammaRayTelAnalysis->SetHisto2DSave(newValue);}

  if( command == Histo2DModeCmd )
    { GammaRayTelAnalysis->SetHisto2DMode(newValue);}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#endif





