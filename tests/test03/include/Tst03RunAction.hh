// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef Tst03RunAction_h
#define Tst03RunAction_h 1

#include "globals.hh"
#include "G4Timer.hh"
#include "G4UserRunAction.hh"

class G4Run;

class Tst03RunAction : public G4UserRunAction
{
  public:

    Tst03RunAction();
    ~Tst03RunAction();

  public:

    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:

    G4Timer* timer;
    G4int runIDcounter;
};

#endif /*Tst03RunAction_h*/
