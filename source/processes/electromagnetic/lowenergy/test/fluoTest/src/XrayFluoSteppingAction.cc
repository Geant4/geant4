//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "XrayFluoSteppingAction.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "G4ios.hh"
#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSteppingAction::XrayFluoSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoSteppingAction::~XrayFluoSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoSteppingAction::UserSteppingAction(const G4Step* aStep)
{
#ifdef G4ANALYSIS_USE
  XrayFluoAnalysisManager* analysis  = XrayFluoAnalysisManager::getInstance();
  analysis->analyseStepping(aStep);
#endif
}	 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
