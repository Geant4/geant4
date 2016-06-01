
#ifndef ExE03RunAction_h
#define ExE03RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "globals.hh"

class G4Run;

class ExE03RunAction : public G4UserRunAction
{
  public:
    ExE03RunAction();
    ~ExE03RunAction();

  public:
    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:
    G4Timer* timer;
    G4int runIDcounter;
};

#endif








