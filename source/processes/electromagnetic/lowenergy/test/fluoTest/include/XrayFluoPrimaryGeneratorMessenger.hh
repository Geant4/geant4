// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: XrayFluoPrimaryGeneratorMessenger.hh,v 1.1 2001-11-27 14:59:32 elena Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoPrimaryGeneratorMessenger_h
#define XrayFluoPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class XrayFluoPrimaryGeneratorAction;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoPrimaryGeneratorMessenger: public G4UImessenger
{
  public:
    XrayFluoPrimaryGeneratorMessenger(XrayFluoPrimaryGeneratorAction*);
   ~XrayFluoPrimaryGeneratorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    XrayFluoPrimaryGeneratorAction* XrayFluoAction; 
    G4UIcmdWithAString*          RndmCmd;
  G4UIcmdWithAString*          RndmVert;
  G4UIcmdWithAString*        spectrum;
  G4UIcmdWithAString*        isoVert;
};

#endif

