//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayFluoAnalysisMessenger.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisMessenger.hh"





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoAnalysisMessenger::XrayFluoAnalysisMessenger(XrayFluoAnalysisManager* analysisManager)
  :xrayFluoAnalysis(analysisManager)
  
{ 
  XrayFluoAnalysisDir = new G4UIdirectory("/analysis/");
  XrayFluoAnalysisDir->SetGuidance("analysis control.");
  
  outputFileCommand = new G4UIcmdWithAString("/analysis/outputFile",this);
  outputFileCommand->SetGuidance("specify the name of the output file (lowercase if paw will be used)");
  outputFileCommand->SetGuidance("command /analysis/update MUST be used after this command");
  outputFileCommand->SetParameterName("choice",true);
  outputFileCommand->SetDefaultValue("xrayfluo.hbk");
  outputFileCommand->AvailableForStates(G4State_Idle);
    
  outputFileType = new G4UIcmdWithAString("/analysis/fileType",this);
  outputFileType->SetGuidance("specify the type of the output file");
  outputFileType->SetGuidance("command /analysis/update MUST be used after this command");
  outputFileType->SetParameterName("choice",true);
  outputFileType->SetDefaultValue("hbook");
  outputFileType->SetCandidates("hbook xml");
  outputFileType->AvailableForStates(G4State_Idle);
  
  persistencyUpdateCommand = new G4UIcmdWithoutParameter("/analysis/update",this);
  outputFileCommand->SetGuidance("Update persistency file");
  outputFileCommand->SetGuidance("This command MUST be used after outputFile or outputFileType commands");
  outputFileCommand->AvailableForStates(G4State_Idle);
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
      if ((xrayFluoAnalysis->GetDeletePersistencyFileFlag())) remove("xrayfluo.hbk");      
      xrayFluoAnalysis->SetOutputFileName(newValue);
    }
  
  else if (command == persistencyUpdateCommand)
    {xrayFluoAnalysis->CreatePersistency();}

  else if(command == outputFileType)
    {
      if (xrayFluoAnalysis->GetDeletePersistencyFileFlag()) remove("xrayfluo.hbk");      
      xrayFluoAnalysis->SetOutputFileType(newValue);
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif













