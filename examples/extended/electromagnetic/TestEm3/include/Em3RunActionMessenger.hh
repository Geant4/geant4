// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3RunActionMessenger.hh,v 1.2 1999-12-15 14:49:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3RunActionMessenger_h
#define Em3RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em3RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3RunActionMessenger: public G4UImessenger
{
  public:
    Em3RunActionMessenger(Em3RunAction*);
   ~Em3RunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em3RunAction*          Em3Run;   
    G4UIcmdWithAString*    SaveCmd;
    G4UIdirectory*         RndmDir;
    G4UIcmdWithAnInteger*  RndmSaveCmd;       
    G4UIcmdWithAString*    RndmReadCmd;
 
};

#endif
