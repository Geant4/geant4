#ifndef BrachyRunAction_h
#define BrachyRunAction_h 1

#include "G4UserRunAction.hh"
#include "G4RunManager.hh"
#include "globals.hh"


class G4Run;
class BrachyAnalysisManager;
class BrachyDetectorConstruction;
class BrachyRunMessenger;

class BrachyFactory;
class BrachyFactoryIr;
class BrachyFactoryI;
class BrachyRunAction : public G4UserRunAction
{
  public:
    BrachyRunAction(G4String& );
   ~BrachyRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run* );
  void SelectEnergy(G4int); 
 private:
  
   G4String      SDname;
   BrachyFactory  *factory; 
   G4VUserDetectorConstruction* pDetector;
   BrachyDetectorConstruction* pDet;
 
  G4int a;
  BrachyRunMessenger* pRun;
 
};

#endif



