//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoEventAction_h
#define XrayFluoEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "XrayFluoAnalysisManager.hh"
#endif


class XrayFluoEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoEventAction : public G4UserEventAction
{
  public:
#ifdef G4ANALYSIS_USE
  XrayFluoEventAction(XrayFluoAnalysisManager* analysisMgr);
#else
    XrayFluoEventAction();
#endif   
 virtual ~XrayFluoEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void SetDrawFlag   (G4String val)  {drawFlag = val;};
    void SetPrintModulo(G4int    val)  {printModulo = val;};
    
  private:
    G4int                       sensorCollID;
    G4int                       sampleCollID;               
    G4String                    drawFlag;
    G4int                       printModulo;                         
    XrayFluoEventActionMessenger*  eventMessenger;
#ifdef G4ANALYSIS_USE
    XrayFluoAnalysisManager* analysisManager;
#endif
};

#endif

    
