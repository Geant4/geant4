// Em6EventAction.hh

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Em6EventAction_h
#define Em6EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class Em6RunAction;
class Em6EventActionMessenger;

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Em6EventAction : public G4UserEventAction
{
  public:

    Em6EventAction(Em6RunAction*);
   ~Em6EventAction();

    void BeginOfEventAction(const G4Event*);
    void   EndOfEventAction(const G4Event*);

    void SetDrawFlag   (G4String val)  {drawFlag    = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};

  private:

    Em6RunAction*             Em6Run;
    G4String                  drawFlag;
    G4int                     printModulo;
    Em6EventActionMessenger*  eventMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


