// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst09RunAction.hh,v 1.4 1999-12-15 14:54:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst09RunAction_h
#define Tst09RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst09RunAction : public G4UserRunAction
{
  public:
    Tst09RunAction();
    virtual ~Tst09RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

