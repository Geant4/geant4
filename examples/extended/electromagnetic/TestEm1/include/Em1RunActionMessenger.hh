// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1RunActionMessenger.hh,v 1.1 1999-10-11 13:07:42 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em1RunActionMessenger_h
#define Em1RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em1RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1RunActionMessenger: public G4UImessenger
{
  public:
    Em1RunActionMessenger(Em1RunAction*);
   ~Em1RunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em1RunAction*          Em1Run;   
    G4UIdirectory*         RndmDir;
    G4UIcmdWithAnInteger*  RndmSaveCmd;    
    G4UIcmdWithAString*    RndmReadCmd;
    
};

#endif
