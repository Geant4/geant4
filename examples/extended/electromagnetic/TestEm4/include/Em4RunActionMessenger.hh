// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em4RunActionMessenger.hh,v 1.1 1999-10-12 11:26:57 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em4RunActionMessenger_h
#define Em4RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em4RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em4RunActionMessenger: public G4UImessenger
{
  public:
    Em4RunActionMessenger(Em4RunAction*);
   ~Em4RunActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em4RunAction*          Em4Run;   
    G4UIcmdWithAString*    SaveCmd;
    G4UIdirectory*         RndmDir;
    G4UIcmdWithAnInteger*  RndmSaveCmd;    
    G4UIcmdWithAString*    RndmReadCmd;
    
};

#endif
