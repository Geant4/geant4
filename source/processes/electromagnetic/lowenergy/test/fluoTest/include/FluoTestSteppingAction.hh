//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestSteppingAction_h
#define FluoTestSteppingAction_h 1

#include "G4UserSteppingAction.hh"


class FluoTestDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestSteppingAction : public G4UserSteppingAction
{
  public:

  FluoTestSteppingAction();   

   ~FluoTestSteppingAction();

    void UserSteppingAction(const G4Step*);
 private:
 

};

#endif
