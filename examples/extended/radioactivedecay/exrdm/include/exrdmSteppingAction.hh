//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef exrdmSteppingAction_h
#define exrdmSteppingAction_h 1

#include "G4UserSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exrdmSteppingAction : public G4UserSteppingAction
{
  public:
    exrdmSteppingAction();
    virtual ~exrdmSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif
