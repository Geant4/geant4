//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestEventAction_h
#define fluoTestEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "fluoTestAnalysisManager.hh"
#endif


class fluoTestEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestEventAction : public G4UserEventAction
{
  public:
#ifdef G4ANALYSIS_USE
  fluoTestEventAction(fluoTestAnalysisManager* analysisMgr);
#else
    fluoTestEventAction();
#endif   
 virtual ~fluoTestEventAction();

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
    fluoTestEventActionMessenger*  eventMessenger;
#ifdef G4ANALYSIS_USE
    fluoTestAnalysisManager* analysisManager;
#endif
};

#endif

    
