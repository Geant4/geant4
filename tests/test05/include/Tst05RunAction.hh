
#ifndef Tst05RunAction_h
#define Tst05RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "globals.hh"

class G4Run;

class Tst05RunAction : public G4UserRunAction
{
  public:
    Tst05RunAction();
    ~Tst05RunAction();

  public:
    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:
    G4Timer* timer;
    G4int runIDcounter;
};

#endif








