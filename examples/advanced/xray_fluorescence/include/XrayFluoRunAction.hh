//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoRunAction_h
#define XrayFluoRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class XrayFluoRunAction : public G4UserRunAction
{
public:
  
#ifdef G4ANALYSIS_USE
  XrayFluoRunAction(XrayFluoAnalysisManager* analysisMgr);
#else 
  XrayFluoRunAction();
#endif 
  ~XrayFluoRunAction();
  
public:
  void BeginOfRunAction(const G4Run*);
  void EndOfRunAction(const G4Run*);
private:
#ifdef G4ANALYSIS_USE
  XrayFluoAnalysisManager* analysisManager;
#endif  
  
};

#endif

