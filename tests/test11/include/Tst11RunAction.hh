// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst11RunAction.hh,v 1.4 1999-12-15 14:54:44 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst11RunAction_h
#define Tst11RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class Tst11RunAction : public G4UserRunAction
{
  public:
    Tst11RunAction();
    virtual ~Tst11RunAction();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

};

#endif

