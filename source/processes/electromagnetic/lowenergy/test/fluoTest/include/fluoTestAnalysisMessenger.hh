//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef fluoTestAnalysisMessenger_h
#define fluoTestAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class fluoTestAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestAnalysisMessenger: public G4UImessenger
{
public:
  fluoTestAnalysisMessenger(fluoTestAnalysisManager* );
  ~fluoTestAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  fluoTestAnalysisManager* fluoTestAnalysis;
  G4UIdirectory*              fluoTestAnalysisDir;
  
  G4UIcmdWithAString*        Histo1DDrawCmd;
  G4UIcmdWithAString*        Histo1DSaveCmd;
 
};
#endif
#endif
