//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef myRunAction_h
#define myRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "myAnalysisManager.hh"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class myRunAction : public G4UserRunAction
{
  public:
   
#ifdef G4ANALYSIS_USE
    myRunAction(myAnalysisManager* analysisMgr);
#else 
   myRunAction();
#endif 
 ~myRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
private:
#ifdef G4ANALYSIS_USE
    myAnalysisManager* analysisManager;
#endif  

};

#endif

