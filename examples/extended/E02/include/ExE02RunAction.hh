
#ifndef ExE02RunAction_h
#define ExE02RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class ExE02RunAction : public G4UserRunAction
{
  public:
    ExE02RunAction();
    ~ExE02RunAction();

  public:
    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

