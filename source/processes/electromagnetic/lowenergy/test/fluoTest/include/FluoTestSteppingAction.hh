//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestSteppingAction_h
#define FluoTestSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#ifdef G4ANALYSIS_USE
#include "FluoTestAnalysisManager.hh"
#endif

class FluoTestDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestSteppingAction : public G4UserSteppingAction
{
  public:
#ifdef G4ANALYSIS_USE    
  //  FluoTestSteppingAction(FluoTestDetectorConstruction*,FluoTestAnalysisManager*);
 FluoTestSteppingAction(FluoTestAnalysisManager*);
#else
  FluoTestSteppingAction();   
#endif
   ~FluoTestSteppingAction();

    void UserSteppingAction(const G4Step*);
 private:
  G4int nElectrons;
  G4int nDepElec;
  G4int nElecCreated;

  //FluoTestDetectorConstruction* detector;
#ifdef G4ANALYSIS_USE   
    FluoTestAnalysisManager* analysisManager;
#endif
};

#endif
