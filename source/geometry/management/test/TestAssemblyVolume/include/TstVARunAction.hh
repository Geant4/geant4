
#ifndef TstVARunAction_h
#define TstVARunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class TstVARunAction : public G4UserRunAction
{
  public:
    TstVARunAction();
    virtual ~TstVARunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

