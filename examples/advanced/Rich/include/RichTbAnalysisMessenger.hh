// Rich advanced example for Geant4
// RichTbAnalysisMessenger.hh for Rich of LHCb
// History:
// Created: Patricia Mendez (Patricia.Mendez@cern.ch)
////////////////////////////////////////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef RichTbAnalysisMessenger_h
#define RichTbAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class RichTbAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class G4UIcmdWithAString;
class RichTbAnalysisMessenger: public G4UImessenger

{
public:
  RichTbAnalysisMessenger(RichTbAnalysisManager* );
  ~RichTbAnalysisMessenger();
  

  
private:

  //pointer to RichTbAnalysisManager
  RichTbAnalysisManager* richAnalysis;
  G4UIdirectory* RichTbAnalysisDir;
  G4UIcmdWithAString* ouputFileCommand;

};
#endif
#endif


