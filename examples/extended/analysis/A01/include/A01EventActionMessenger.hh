// $Id: A01EventActionMessenger.hh,v 1.1 2002-11-13 07:18:10 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#ifndef A01EventActionMessenger_h
#define A01EventActionMessenger_h 1

class A01EventAction;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

class A01EventActionMessenger: public G4UImessenger
{
  public:
    A01EventActionMessenger(A01EventAction* mpga);
    ~A01EventActionMessenger();

  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    A01EventAction* target;

  private: //commands
    G4UIcmdWithAnInteger*  verboseCmd;

};

#endif


