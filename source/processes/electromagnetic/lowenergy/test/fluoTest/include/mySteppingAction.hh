//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef mySteppingAction_h
#define mySteppingAction_h 1

#include "G4UserSteppingAction.hh"
#ifdef G4ANALYSIS_USE
#include "myAnalysisManager.hh"
#endif

class myDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class mySteppingAction : public G4UserSteppingAction
{
  public:
#ifdef G4ANALYSIS_USE    
  mySteppingAction(myDetectorConstruction*,myAnalysisManager*);
#else
  mySteppingAction();   
#endif
   ~mySteppingAction();

    void UserSteppingAction(const G4Step*);
 private:
    myDetectorConstruction* detector;
#ifdef G4ANALYSIS_USE   
    myAnalysisManager* analysisManager;
#endif
};

#endif
