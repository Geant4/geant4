//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestEventAction_h
#define FluoTestEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
class FluoTestAnalysisManager;

#endif


class FluoTestEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestEventAction : public G4UserEventAction
{
  public:

    FluoTestEventAction();
   
   virtual ~FluoTestEventAction();

    public:
      virtual void   BeginOfEventAction(const G4Event*);
      virtual void   EndOfEventAction(const G4Event*);
    
     void SetDrawFlag   (G4String val)  {drawFlag = val;};
      void SetPrintModulo(G4int    val)  {printModulo = val;};
    
  private:
  //  G4int                       HPGeCollID;
   G4String                    drawFlag;
    G4int                       printModulo;                         
    FluoTestEventActionMessenger*  eventMessenger;

};

#endif

    
