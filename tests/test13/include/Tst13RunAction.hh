// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst13RunAction.hh,v 1.1 1999-01-08 16:35:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst13RunAction_h
#define Tst13RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst13RunAction : public G4UserRunAction
{
  public:
    Tst13RunAction();
    ~Tst13RunAction();

  public:
    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

