
#ifndef Tst01RunAction_h
#define Tst01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst01RunAction : public G4UserRunAction
{
  public:
    Tst01RunAction();
    ~Tst01RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

