// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst01RunAction.hh,v 1.1 2001-02-08 08:41:45 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst01RunAction_h
#define Tst01RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Timer.hh"

class G4Run;

class Tst01RunAction : public G4UserRunAction
{
  public:
    Tst01RunAction();
    ~Tst01RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);
  
  private:
    G4Timer  ftimer;
};

#endif

