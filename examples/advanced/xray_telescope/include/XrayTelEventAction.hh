//  XrayTelEventAction.hh

#ifndef XrayTelEventAction_h
#define XrayTelEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class XrayTelEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayTelEventAction : public G4UserEventAction
{
  public:
    XrayTelEventAction();
   ~XrayTelEventAction();

  public:
    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    
  private:
    G4String drawFlag;                         // control the drawing of event
    XrayTelEventActionMessenger*  eventMessenger;
};

#endif

    




