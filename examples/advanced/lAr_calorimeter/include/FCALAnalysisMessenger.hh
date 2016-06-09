// Author: Patricia Mendez (patricia.mendez@cern.ch)
//
// History:
// -----------
//  14 Feb 2003  Patricia Mendez   Created
//
// -------------------------------------------------------------------




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef FCALAnalysisMessenger_h
#define FCALAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class FCALAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class G4UIcmdWithAString;
class FCALAnalysisMessenger: public G4UImessenger

{
public:
  FCALAnalysisMessenger(FCALAnalysisManager* );
  ~FCALAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  //pointer to FCALAnalysisManager
  FCALAnalysisManager* xrayFluoAnalysis;
  G4UIdirectory* FCALAnalysisDir;
  G4UIcmdWithAString* ouputFileCommand;

};
#endif
#endif


