// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StackingMessenger.hh,v 1.3 1999-12-15 14:49:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4StackingMessenger_h
#define G4StackingMessenger_h 1

#include "G4UImessenger.hh"
class G4StackManager;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger;

// class description:
//
//  This is a concrete class of G4UImessenger which handles the commands
// for G4StackManager. It has the following commands:
//   /event/stack/
//   /event/stack/status
//   /event/stack/clear

class G4StackingMessenger: public G4UImessenger
{
  public:
    G4StackingMessenger(G4StackManager* fCont);
    ~G4StackingMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    G4StackManager * fContainer;
    G4UIdirectory* stackDir;
    G4UIcmdWithoutParameter* statusCmd;
    G4UIcmdWithAnInteger* clearCmd;
};

#endif

