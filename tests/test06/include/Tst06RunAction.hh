// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst06RunAction.hh,v 1.1 1999-01-08 16:35:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst06RunAction_h
#define Tst06RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst06RunAction : public G4UserRunAction
{
  public:
    Tst06RunAction();
    ~Tst06RunAction();

  public:
    void BeginOfRunAction(G4Run* aRun);
    void EndOfRunAction(G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

