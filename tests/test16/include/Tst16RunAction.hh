// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst16RunAction.hh,v 1.2 1999-12-15 14:54:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst16RunAction_h
#define Tst16RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst16RunAction : public G4UserRunAction
{
  public:
    Tst16RunAction();
    ~Tst16RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

};

#endif

