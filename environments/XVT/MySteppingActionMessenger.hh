// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingActionMessenger.hh,v 2.1 1998/07/12 02:37:20 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef MySteppingActionMessenger_h
#define MySteppingActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"

class MySteppingAction;

class MySteppingActionMessenger: public G4UImessenger
{
  public:
    MySteppingActionMessenger(MySteppingAction * mySA);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    MySteppingAction * mySteppingAction;
};

#endif

