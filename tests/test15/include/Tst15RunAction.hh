// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst15RunAction.hh,v 1.2 1999-12-15 14:54:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst15RunAction_h
#define Tst15RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst15RunAction : public G4UserRunAction
{
  public:
    Tst15RunAction();
    ~Tst15RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

};

#endif

