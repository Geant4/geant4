//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoSteppingAction_h
#define XrayFluoSteppingAction_h 1
#include "globals.hh"
#include "G4UserSteppingAction.hh"

class XrayFluoDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoSteppingAction : public G4UserSteppingAction
{
  public:

  XrayFluoSteppingAction(); 
   ~XrayFluoSteppingAction();

    void UserSteppingAction(const G4Step*);
private:


};

#endif
