// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4RunAction.hh,v 1.3 1999/12/05 17:50:26 gcosmo Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef G3toG4RunAction_h
#define G3toG4RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class G3toG4RunAction : public G4UserRunAction
{
  public:
    G3toG4RunAction();
    ~G3toG4RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif

