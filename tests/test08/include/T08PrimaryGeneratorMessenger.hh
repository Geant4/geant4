// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08PrimaryGeneratorMessenger.hh,v 1.1 1999-01-08 16:35:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef T08PrimaryGeneratorMessenger_h
#define T08PrimaryGeneratorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class T08PrimaryGeneratorAction;
class G4UIcmdWithAString;

class T08PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    T08PrimaryGeneratorMessenger(T08PrimaryGeneratorAction* myGun);
    ~T08PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand * command,G4String newValues);
    
  private:
    T08PrimaryGeneratorAction* myAction;
 
    G4UIcmdWithAString*       RndmCmd;
};

#endif

