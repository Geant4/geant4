
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef GammaRayTelEventAction_h
#define GammaRayTelEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//class GammaRayTelEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelEventAction : public G4UserEventAction
{
  public:
    GammaRayTelEventAction();
    virtual ~GammaRayTelEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};

  private:
    G4int       trackerCollID;                
    G4String    drawFlag;
    G4int       printModulo;                         
};

#endif

    



