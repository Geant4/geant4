//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestRunAction_h
#define FluoTestRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class FluoTestRunAction : public G4UserRunAction
{
  public:
   
 
   FluoTestRunAction();

 ~FluoTestRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
private:
  

};

#endif

