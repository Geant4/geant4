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













