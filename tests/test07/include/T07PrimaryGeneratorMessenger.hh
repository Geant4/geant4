// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07PrimaryGeneratorMessenger.hh,v 1.1 1999-01-08 16:35:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef T07PrimaryGeneratorMessenger_h
#define T07PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class T07PrimaryGeneratorAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class T07PrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    T07PrimaryGeneratorMessenger(T07PrimaryGeneratorAction*);
   ~T07PrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    T07PrimaryGeneratorAction* T07Action; 
    G4UIcmdWithAString*          RndmCmd;
};

#endif

