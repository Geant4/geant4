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
// Author: Patricia Mendez (patricia.mendez@cern.ch)
//
// History:
// -----------
// 14 Feb 2003 Patricia Mendez     Created
//--------------------------------------------------------------------------

#ifdef G4ANALYSIS_USE
#include "FCALAnalysisMessenger.hh"

#include "FCALAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisMessenger::FCALAnalysisMessenger(FCALAnalysisManager* analysisManager)
  :xrayFluoAnalysis(analysisManager)

{ 
  FCALAnalysisDir = new G4UIdirectory("/analysis/");
  FCALAnalysisDir->SetGuidance("analysis control.");

  ouputFileCommand = new G4UIcmdWithAString("/analysis/outputFile",this);
  ouputFileCommand->SetGuidance("specify the name of the output file (lowercase for Hbook)");
  ouputFileCommand->SetGuidance("default: xrayfluo.hbk");
  ouputFileCommand->SetParameterName("choice",true);
  ouputFileCommand->SetDefaultValue("fcal.his");
  ouputFileCommand->AvailableForStates(G4State_Idle);
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALAnalysisMessenger::~FCALAnalysisMessenger()
{
  
  delete FCALAnalysisDir; 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALAnalysisMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

if(command == ouputFileCommand)
    {xrayFluoAnalysis->SetOutputFileName(newValue);}

   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif













