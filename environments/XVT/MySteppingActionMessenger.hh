// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingActionMessenger.hh,v 1.2 1999-12-15 14:48:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

