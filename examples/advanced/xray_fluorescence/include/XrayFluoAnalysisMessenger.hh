//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
#ifndef XrayFluoAnalysisMessenger_h
#define XrayFluoAnalysisMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class XrayFluoAnalysisManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoAnalysisMessenger: public G4UImessenger
{
public:
  XrayFluoAnalysisMessenger(XrayFluoAnalysisManager* );
  ~XrayFluoAnalysisMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  XrayFluoAnalysisManager* XrayFluoAnalysis;
  G4UIdirectory*              XrayFluoAnalysisDir;
  
  G4UIcmdWithAString*        Histo1DDrawCmd;
  G4UIcmdWithAString*        Histo1DSaveCmd;
 
};
#endif
#endif

