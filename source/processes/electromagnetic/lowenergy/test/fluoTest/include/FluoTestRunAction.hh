//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestRunAction_h
#define FluoTestRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "FluoTestAnalysisManager.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class FluoTestRunAction : public G4UserRunAction
{
  public:
   
#ifdef G4ANALYSIS_USE
    FluoTestRunAction(FluoTestAnalysisManager* analysisMgr);
#else 
   FluoTestRunAction();
#endif 
 ~FluoTestRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
private:
#ifdef G4ANALYSIS_USE
    FluoTestAnalysisManager* analysisManager;
#endif  

};

#endif

