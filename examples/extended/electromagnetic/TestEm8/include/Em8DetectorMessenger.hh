// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8DetectorMessenger.hh,v 1.1 2000-01-07 14:50:21 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em8DetectorMessenger_h
#define Em8DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em8DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em8DetectorMessenger: public G4UImessenger
{
  public:
    Em8DetectorMessenger(Em8DetectorConstruction* );
   ~Em8DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em8DetectorConstruction*   Em8Detector;
    
    G4UIdirectory*             Em8detDir;

    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsRadCmd;

    G4UIcmdWithADoubleAndUnit* AbsZposCmd;

    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldZCmd;
    G4UIcmdWithADoubleAndUnit* WorldRCmd;

    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

};

#endif

