
#ifndef PersEx01RunAction_h
#define PersEx01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class PersEx01RunAction : public G4UserRunAction
{
  public:
    PersEx01RunAction();
    virtual ~PersEx01RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
};

#endif

