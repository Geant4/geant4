//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestRunAction_h
#define fluoTestRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisManager.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class fluoTestRunAction : public G4UserRunAction
{
  public:
   
#ifdef G4ANALYSIS_USE
    fluoTestRunAction(fluoTestAnalysisManager* analysisMgr);
#else 
   fluoTestRunAction();
#endif 
 ~fluoTestRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
private:
#ifdef G4ANALYSIS_USE
    fluoTestAnalysisManager* analysisManager;
#endif  

};

#endif

