
#ifndef exGPSEventAction_h
#define exGPSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class exGPSEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exGPSEventAction : public G4UserEventAction
{
  public:
    exGPSEventAction() ;
    virtual ~exGPSEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    
  private:
    G4String                    drawFlag;
    G4int                       printModulo;                         
    exGPSEventActionMessenger*  eventMessenger;
};

#endif

    
