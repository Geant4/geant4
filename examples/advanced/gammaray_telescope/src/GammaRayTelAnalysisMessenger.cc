//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: GammaRayTelAnalysisMessenger.cc,v 1.9 2006-06-29 15:56:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
//
//      ------------ GammaRayTelAnalysisMessenger ------
//           by R.Giannitrapani, F.Longo & G.Santin (03 dic 2000)
// 20.11.01 G.Santin: modified according to the new GammaRayTelAnalysis.cc
// ************************************************************
#ifdef G4ANALYSIS_USE

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

#include "GammaRayTelAnalysisMessenger.hh"
#include "GammaRayTelAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GammaRayTelAnalysisMessenger::GammaRayTelAnalysisMessenger(GammaRayTelAnalysis* analysis)
  :gammaRayTelAnalysis(analysis)

{ 
  gammaRayTelAnalysisDir = new G4UIdirectory("/analysis/");
  gammaRayTelAnalysisDir->SetGuidance("GammaRayTel analysis control.");
  
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
  Histo1DDrawCmd->AvailableForStates(G4State_Idle);

  Histo1DSaveCmd = new G4UIcmdWithAString("/analysis/histo1dSave",this);
  Histo1DSaveCmd->SetGuidance("Enable the saving of the 1d histograms every run.");
  Histo1DSaveCmd->SetGuidance("Choice: disable, enable(default)");
  Histo1DSaveCmd->SetParameterName("choice",true);
  Histo1DSaveCmd->SetDefaultValue("enable");
  Histo1DSaveCmd->SetCandidates("disable enable");
  Histo1DSaveCmd->AvailableForStates(G4State_Idle);

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
  Histo2DDrawCmd->AvailableForStates(G4State_Idle);

  Histo2DSaveCmd = new G4UIcmdWithAString("/analysis/histo2dSave",this);
  Histo2DSaveCmd->SetGuidance("Enable the saving of the 2d histograms every run.");
  Histo2DSaveCmd->SetGuidance("Choice: disable, enable(default)");
  Histo2DSaveCmd->SetParameterName("choice",true);
  Histo2DSaveCmd->SetDefaultValue("enable");
  Histo2DSaveCmd->SetCandidates("disable enable");
  Histo2DSaveCmd->AvailableForStates(G4State_Idle);

  Histo2DModeCmd = new G4UIcmdWithAString("/analysis/histo2dMode",this);
  Histo2DModeCmd->SetGuidance("Select the mode for the 2d histograms.");
  Histo2DModeCmd->SetGuidance("Choice: position, strip(default)");
  Histo2DModeCmd->SetGuidance("position -> the histo is filled with true positions in mm");
  Histo2DModeCmd->SetGuidance("strip -> the histo is filled with the number of the strip and the plane");
  Histo2DModeCmd->SetParameterName("choice",true);
  Histo2DModeCmd->SetDefaultValue("strip");
  Histo2DModeCmd->SetCandidates("position strip");
  Histo2DModeCmd->AvailableForStates(G4State_Idle);
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
    { gammaRayTelAnalysis->SetHisto1DDraw(newValue);}

  if( command == Histo1DSaveCmd )
    { gammaRayTelAnalysis->SetHisto1DSave(newValue);}

  // 2D Histograms

  if( command == Histo2DDrawCmd )
    { gammaRayTelAnalysis->SetHisto2DDraw(newValue);}

  if( command == Histo2DSaveCmd )
    { gammaRayTelAnalysis->SetHisto2DSave(newValue);}

  if( command == Histo2DModeCmd )
    { gammaRayTelAnalysis->SetHisto2DMode(newValue);}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#endif





