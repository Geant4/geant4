
#ifndef ThyroidRunAction_h
#define ThyroidRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"


class G4Run;
class ThyroidAnalysisManager;
class ThyroidRunAction : public G4UserRunAction
{
  public:
    ThyroidRunAction( );
   ~ThyroidRunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run* );
  
};

#endif


