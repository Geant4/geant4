//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoEventAction_h
#define XrayFluoEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class XrayFluoRunAction;
class XrayFluoEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


class XrayFluoEventAction : public G4UserEventAction
{
public:
  
  XrayFluoEventAction();
  
  virtual ~XrayFluoEventAction();
  
public:
  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  
  void SetDrawFlag   (G4String val)  {drawFlag = val;};
  void SetPrintModulo(G4int    val)  {printModulo = val;};
  
private:
  
  G4String                    drawFlag;
  G4int                       HPGeCollID; 
  XrayFluoEventActionMessenger*  eventMessenger;
  G4int                       printModulo;                         
  
  G4double RandomCut(G4double);
  G4double ResponseFunction(G4double);
  
  XrayFluoRunAction* runManager;
};

#endif

    
