//  XrayTelRunAction

#ifndef XrayTelRunAction_h
#define XrayTelRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4Run;

class XrayTelRunAction : public G4UserRunAction
{
  public:
    XrayTelRunAction();
   ~XrayTelRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

};

#endif

