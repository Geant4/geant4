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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoAnalysisMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisMessenger::XrayFluoAnalysisMessenger(XrayFluoAnalysisManager* analysisManager)
  :xrayFluoAnalysis(analysisManager)
  
{ 
  XrayFluoAnalysisDir = new G4UIdirectory("/analysis/");
  XrayFluoAnalysisDir->SetGuidance("analysis control.");
  
  outputFileCommand = new G4UIcmdWithAString("/analysis/outputFile",this);
  outputFileCommand->SetGuidance("specify the name of the output file");
  outputFileCommand->SetParameterName("choice",true);
  outputFileCommand->SetDefaultValue("xrayfluo.root");
  outputFileCommand->AvailableForStates(G4State_Idle);
    

  physicFlagCmd = new G4UIcmdWithABool("/analysis/setPhysicProduction",this);
  physicFlagCmd->SetGuidance("Select if data stored in the Phase-Space must contain physical data or particles exiting the sample");
  physicFlagCmd->SetGuidance("To be used before and togheter with /gun/loadGunData");
  physicFlagCmd->SetParameterName("Physyc Flag",true);
  physicFlagCmd->SetDefaultValue(false);
  physicFlagCmd->AvailableForStates(G4State_Idle);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisMessenger::~XrayFluoAnalysisMessenger()
{
  
  delete XrayFluoAnalysisDir; 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if(command == outputFileCommand)
    {      
      xrayFluoAnalysis->SetOutputFileName(newValue);
    }
  
  if( command == physicFlagCmd )
    { 
      G4bool newPhysFlag = physicFlagCmd->GetNewBoolValue(newValue);
      xrayFluoAnalysis->SetPhysicFlag(newPhysFlag);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....















