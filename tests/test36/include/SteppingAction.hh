#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(DetectorConstruction*, RunAction*);
   ~SteppingAction();

    void UserSteppingAction(const G4Step*);
    
    
  private:
    DetectorConstruction* detector;
    RunAction*          runaction;      

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
