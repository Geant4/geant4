//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestSteppingAction_h
#define fluoTestSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisManager.hh"
#endif

class fluoTestDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestSteppingAction : public G4UserSteppingAction
{
  public:
#ifdef G4ANALYSIS_USE    
  fluoTestSteppingAction(fluoTestDetectorConstruction*,fluoTestAnalysisManager*);
#else
  fluoTestSteppingAction();   
#endif
   ~fluoTestSteppingAction();

    void UserSteppingAction(const G4Step*);
 private:
    fluoTestDetectorConstruction* detector;
#ifdef G4ANALYSIS_USE   
    fluoTestAnalysisManager* analysisManager;
#endif
};

#endif
