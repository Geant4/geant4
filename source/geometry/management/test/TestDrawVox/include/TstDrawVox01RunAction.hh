
#ifndef TstDrawVox01RunAction_h
#define TstDrawVox01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class TstDrawVox01RunAction : public G4UserRunAction
{
  public:
    TstDrawVox01RunAction();
    virtual ~TstDrawVox01RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

