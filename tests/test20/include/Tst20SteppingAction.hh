// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef Tst20SteppingAction_H
#define Tst20SteppingAction_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include "Tst20ProcCallSA.hh"

class Tst20SteppingAction : public G4UserSteppingAction
{
  public:
    Tst20SteppingAction();
    virtual ~Tst20SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
    
  private:
    Tst20ProcCallSA  theProcCallSA;
};

#endif

