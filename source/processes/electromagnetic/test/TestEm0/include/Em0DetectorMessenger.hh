// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0DetectorMessenger.hh,v 1.1 1999-01-08 16:32:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em0DetectorMessenger_h
#define Em0DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em0DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0DetectorMessenger: public G4UImessenger
{
  public:
  
    Em0DetectorMessenger(Em0DetectorConstruction* );
   ~Em0DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
    Em0DetectorConstruction*   Em0Detector;
    
    G4UIdirectory*             Em0detDir;
    G4UIcmdWithAString*        MaterCmd;
    G4UIcmdWithADoubleAndUnit* SizeCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif

