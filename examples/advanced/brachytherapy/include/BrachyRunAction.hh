#ifndef BrachyRunAction_h
#define BrachyRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


class G4Run;
class BrachyAnalysisManager;
class BrachyRunAction : public G4UserRunAction
{
  public:
    BrachyRunAction(G4String& );
   ~BrachyRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run* );
  private:
 
   G4String      SDname;
  
 
};

#endif


