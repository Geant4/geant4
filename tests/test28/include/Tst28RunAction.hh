
#ifndef Tst28RunAction_h
#define Tst28RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst28RunAction : public G4UserRunAction
{
  public:
    Tst28RunAction();
    virtual ~Tst28RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

