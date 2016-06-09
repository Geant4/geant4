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
// RichTbAnalysisMessenger.cc for Rich of LHCb
// Advanced example of Geant4
// Author: Patricia Mendez (patricia.mendez@cern.ch)
//
// History:
// -----------
// Patricia Mendez     Created
//--------------------------------------------------------------------------

#ifdef G4ANALYSIS_USE
#include "RichTbAnalysisMessenger.hh"

#include "RichTbAnalysisManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RichTbAnalysisMessenger::RichTbAnalysisMessenger(RichTbAnalysisManager* analysisManager)
  :richAnalysis(analysisManager)

{ 
  RichTbAnalysisDir = new G4UIdirectory("/analysis/");
  RichTbAnalysisDir->SetGuidance("analysis control.");

  ouputFileCommand = new G4UIcmdWithAString("/analysis/outputFile",this);
  ouputFileCommand->SetGuidance("specify the name of the output file (lowercase for Hbook)");
  ouputFileCommand->SetGuidance("default: rich.hbk");
  ouputFileCommand->SetParameterName("choice",true);
  ouputFileCommand->SetDefaultValue("rich.his");
  ouputFileCommand->AvailableForStates(G4State_Idle);
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RichTbAnalysisMessenger::~RichTbAnalysisMessenger()
{
  
  delete RichTbAnalysisDir; 
 
}

#endif













