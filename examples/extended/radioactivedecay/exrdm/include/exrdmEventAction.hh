
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef exrdmEventAction_h
#define exrdmEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class exrdmEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exrdmEventAction : public G4UserEventAction
{
  public:
    exrdmEventAction();
   ~exrdmEventAction();

  public:
    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    
  private:
    G4String drawFlag;                         // control the drawing of event
    exrdmEventActionMessenger*  eventMessenger;
};

#endif

    




