// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20RunAction.hh,v 1.2 2001-05-25 12:50:06 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst20RunAction_h
#define Tst20RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst20RunAction : public G4UserRunAction
{
  public:
    Tst20RunAction();
    virtual ~Tst20RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

