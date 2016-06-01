// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyEventActionMessenger.hh,v 2.1 1998/07/12 02:37:17 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef MyEventActionMessenger_h
#define MyEventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"

class MyEventAction;

class MyEventActionMessenger: public G4UImessenger
{
  public:
    MyEventActionMessenger(MyEventAction * myEA);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    MyEventAction * myEventAction;
};

#endif

