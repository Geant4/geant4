// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4RunAction.hh,v 1.1 2000-07-24 11:23:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

