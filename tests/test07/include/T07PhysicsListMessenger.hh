// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07PhysicsListMessenger.hh,v 1.2 1999-12-15 14:54:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef T07PhysicsListMessenger_h
#define T07PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class T07PhysicsList;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class T07PhysicsListMessenger: public G4UImessenger
{
  public:
    T07PhysicsListMessenger(T07PhysicsList*);
   ~T07PhysicsListMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    T07PhysicsList*        T07List;
    G4UIcmdWithoutParameter* ProcessCmd;
};

#endif

