#ifndef exGPSAnalysisMessenger_h
#define exGPSAnalysisMessenger_h 1

#ifdef G4ANALYSIS_USE

#include "globals.hh"
#include "G4UImessenger.hh"

class exGPSAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSAnalysisMessenger: public G4UImessenger
{
public:
  exGPSAnalysisMessenger(exGPSAnalysisManager* );
  ~exGPSAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  exGPSAnalysisManager* exGPSAnalysis;
  G4UIdirectory*              exGPSAnalysisDir;
  
  G4UIcmdWithAString*               FileNameCmd;
  G4UIcmdWithAString*               FileTypeCmd;
  G4UIcmdWithADoubleAndUnit*        MaxEngCmd;
  G4UIcmdWithADoubleAndUnit*        MinEngCmd;
  G4UIcmdWithADoubleAndUnit*        MaxPosCmd;
  G4UIcmdWithADoubleAndUnit*        MinPosCmd;
};

#endif // G4ANALYSIS_USE

#endif










