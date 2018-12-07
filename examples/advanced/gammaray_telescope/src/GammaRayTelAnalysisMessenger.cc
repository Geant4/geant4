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
  delete Histo2DModeCmd; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void GammaRayTelAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == Histo2DModeCmd )
    { gammaRayTelAnalysis->SetHisto2DMode(newValue);}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





