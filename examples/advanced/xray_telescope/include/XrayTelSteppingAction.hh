//  XrayTelSteppingAction.hh

#ifndef XrayTelSteppingAction_h
#define XrayTelSteppingAction_h 1

#include "G4UserSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelSteppingAction : public G4UserSteppingAction
{
  public:
    XrayTelSteppingAction();
    virtual ~XrayTelSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
};

#endif
