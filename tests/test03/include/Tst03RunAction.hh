// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef Tst03RunAction_h
#define Tst03RunAction_h 1

#include "G4Timer.hh"
#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Run;

class Tst03RunAction : public G4UserRunAction
{
  public:

    Tst03RunAction();
    virtual ~Tst03RunAction();

  public:

    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:

    G4Timer* timer;
};

#endif /*Tst03RunAction_h*/
