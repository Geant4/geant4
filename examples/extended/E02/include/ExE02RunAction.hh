
#ifndef ExE02RunAction_h
#define ExE02RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class ExE02RunAction : public G4UserRunAction
{
  public:
    ExE02RunAction();
    virtual ~ExE02RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
};

#endif

