//  XrayTelSteppingAction.hh

#ifndef XrayTelSteppingAction_h
#define XrayTelSteppingAction_h 1

class XrayTelHistogram;

#include "G4UserSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelSteppingAction : public G4UserSteppingAction
{
  public:
    XrayTelSteppingAction(XrayTelHistogram* histoMgr);
    virtual ~XrayTelSteppingAction();

    virtual void UserSteppingAction(const G4Step*);
  
  private:
  
    XrayTelHistogram* histoManager;

};

#endif
