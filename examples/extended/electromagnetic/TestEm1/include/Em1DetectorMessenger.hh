// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1DetectorMessenger.hh,v 1.2 1999-12-15 14:48:53 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em1DetectorMessenger_h
#define Em1DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em1DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1DetectorMessenger: public G4UImessenger
{
  public:
  
    Em1DetectorMessenger(Em1DetectorConstruction* );
   ~Em1DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
    Em1DetectorConstruction*   Em1Detector;
    
    G4UIdirectory*             Em1detDir;
    G4UIcmdWithAString*        MaterCmd;
    G4UIcmdWithADoubleAndUnit* SizeCmd;
    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif

