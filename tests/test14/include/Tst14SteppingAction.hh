// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef Tst14SteppingAction_H
#define Tst14SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "Tst14ProcCallSA.hh"

class Tst14SteppingAction : public G4UserSteppingAction
{
  public:
    Tst14SteppingAction();
    virtual ~Tst14SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
    
  private:
    Tst14ProcCallSA  theProcCallSA;
};

#endif

