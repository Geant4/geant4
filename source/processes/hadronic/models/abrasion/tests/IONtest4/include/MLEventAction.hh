#ifndef MLEventAction_h
#define MLEventAction_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "G4UserEventAction.hh"
#include "globals.hh"

#include "MLAnalysisManager.hh"
#include "MLGeometryConstruction.hh"
#include "MLEventActionMessenger.hh"
#include "MLFluenceAnalyser.hh"
#include "MLNIELAnalyser.hh"
#include "MLDoseAnalyser.hh"
#include "MLPHSAnalyser.hh"
////////////////////////////////////////////////////////////////////////////////
//
class MLEventAction : public G4UserEventAction
{
public:
  MLEventAction (MLGeometryConstruction*) ;
  virtual ~MLEventAction();

public:
  virtual void   BeginOfEventAction (const G4Event*);
  virtual void   EndOfEventAction (const G4Event*);

  void SetRandomSeed (G4int);
  void SetCPUTime (G4double dval) {cputime = dval;};
  void SetDrawFlag (G4String sval) {drawFlag = sval;};
  void SetPrintModulo (G4int ival) {printModulo = ival;};

private:
  G4int                    edetectorID;
  MLGeometryConstruction  *geometry;
  MLAnalysisManager       *analysisManager;
  MLEventActionMessenger  *eventMessenger;
  MLFluenceAnalyser       *fluenceAnalyser;
  MLNIELAnalyser          *NIELAnalyser;
  MLDoseAnalyser          *doseAnalyser;
  MLPHSAnalyser           *PHSAnalyser;

  G4String    drawFlag;
  G4int       printModulo;
  G4double    cputime;
};
////////////////////////////////////////////////////////////////////////////////
#endif
