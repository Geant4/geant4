
#ifndef MyRunAction_h
#define MyRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class MyRunAction : public G4UserRunAction
{
  public:
    MyRunAction();
    ~MyRunAction();

  public:
    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

