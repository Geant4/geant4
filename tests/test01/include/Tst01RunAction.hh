
#ifndef Tst01RunAction_h
#define Tst01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst01RunAction : public G4UserRunAction
{
  public:
    Tst01RunAction();
    virtual ~Tst01RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

