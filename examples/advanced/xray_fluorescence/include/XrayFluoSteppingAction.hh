//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoSteppingAction_h
#define  XrayFluoSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif

class  XrayFluoDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class  XrayFluoSteppingAction : public G4UserSteppingAction
{
public:
#ifdef G4ANALYSIS_USE    
  XrayFluoSteppingAction( XrayFluoDetectorConstruction*,XrayFluoAnalysisManager*);
#else
  XrayFluoSteppingAction();   
#endif
  ~ XrayFluoSteppingAction();
  
  void UserSteppingAction(const G4Step*);
private:
  XrayFluoDetectorConstruction* detector;
#ifdef G4ANALYSIS_USE   
  XrayFluoAnalysisManager* analysisManager;
#endif
};

#endif
