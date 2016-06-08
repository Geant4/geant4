
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
    virtual ~ExE03RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
    G4Timer* timer;
};

#endif








