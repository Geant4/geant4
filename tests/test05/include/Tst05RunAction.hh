// $Id: Tst05RunAction.hh,v 1.5 2000-02-25 16:56:41 gcosmo Exp $
// ------------------------------------------------------------

#ifndef Tst05RunAction_h
#define Tst05RunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;

class Tst05RunAction : public G4UserRunAction
{
  public:
    Tst05RunAction();
    virtual ~Tst05RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
    G4Timer* timer;
};

#endif
