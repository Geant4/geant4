// Em6SteppingAction.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6SteppingAction_h
#define Em6SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class Em6DetectorConstruction;
class Em6RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6SteppingAction : public G4UserSteppingAction
{
  public:
    Em6SteppingAction(Em6DetectorConstruction*,Em6RunAction*);
   ~Em6SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
  Em6DetectorConstruction* Em6Det;
  Em6RunAction*            Em6Run;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
