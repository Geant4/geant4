//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef myEventAction_h
#define myEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "myAnalysisManager.hh"
#endif


class myEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class myEventAction : public G4UserEventAction
{
  public:
#ifdef G4ANALYSIS_USE
  myEventAction(myAnalysisManager* analysisMgr);
#else
    myEventAction();
#endif   
 virtual ~myEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    
  private:
    G4int                       siCollID;
    G4int                       sampleCollID;
    G4int                       hpgeCollID;
    G4String                    drawFlag;
    G4int                       printModulo;                         
    myEventActionMessenger*  eventMessenger;
#ifdef G4ANALYSIS_USE
    myAnalysisManager* analysisManager;
#endif
};

#endif

    
