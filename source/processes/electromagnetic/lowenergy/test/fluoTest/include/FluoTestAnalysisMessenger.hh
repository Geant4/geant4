//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef FluoTestAnalysisMessenger_h
#define FluoTestAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class FluoTestAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestAnalysisMessenger: public G4UImessenger
{
public:
  FluoTestAnalysisMessenger(FluoTestAnalysisManager* );
  ~FluoTestAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  FluoTestAnalysisManager* FluoTestAnalysis;
  G4UIdirectory*              FluoTestAnalysisDir;
  
  G4UIcmdWithAString*        Histo1DDrawCmd;
  //G4UIcmdWithAString*        Histo2DDrawCmd;
  G4UIcmdWithAString*        Histo1DSaveCmd;
  //G4UIcmdWithAString*        Histo2DSaveCmd;
  //G4UIcmdWithAString*        Histo2DModeCmd;
};
#endif
#endif
