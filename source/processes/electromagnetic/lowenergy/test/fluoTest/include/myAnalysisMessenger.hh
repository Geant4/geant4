//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef myAnalysisMessenger_h
#define myAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class myAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myAnalysisMessenger: public G4UImessenger
{
public:
  myAnalysisMessenger(myAnalysisManager* );
  ~myAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  myAnalysisManager* myAnalysis;
  G4UIdirectory*              myAnalysisDir;
  
  G4UIcmdWithAString*        Histo1DDrawCmd;
  //G4UIcmdWithAString*        Histo2DDrawCmd;
  G4UIcmdWithAString*        Histo1DSaveCmd;
  //G4UIcmdWithAString*        Histo2DSaveCmd;
  //G4UIcmdWithAString*        Histo2DModeCmd;
};
#endif
#endif
