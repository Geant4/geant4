// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst02RunAction.hh,v 1.2 1999-04-16 09:41:23 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst02RunAction_h
#define Tst02RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst02RunAction : public G4UserRunAction
{
  public:
    Tst02RunAction();
    ~Tst02RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

