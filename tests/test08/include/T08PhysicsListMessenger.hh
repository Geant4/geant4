// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08PhysicsListMessenger.hh,v 1.1 1999-01-08 16:35:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef T08PhysicsListMessenger_h
#define T08PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class T08PhysicsList;
class G4UIcmdWithoutParameter;

class T08PhysicsListMessenger: public G4UImessenger
{
  public:
    T08PhysicsListMessenger(T08PhysicsList* List);
    ~T08PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand* command,G4String newValues);
    
  private:
    T08PhysicsList* myList;
    G4UIcmdWithoutParameter* ProcessCmd;
};

#endif

