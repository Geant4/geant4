#ifndef exGPSRunAction_h
#define exGPSRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class exGPSRunAction : public G4UserRunAction
{
  public:
  exGPSRunAction();
  ~exGPSRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
  private:
};

#endif



