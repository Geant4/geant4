// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// Detector mesenger is defined.
// Class Description - end
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em6DetectorMessenger_h
#define Em6DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Em6DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6DetectorMessenger: public G4UImessenger
{
public: // Without description

    Em6DetectorMessenger(Em6DetectorConstruction* );
   ~Em6DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Em6DetectorConstruction*   Em6Detector;
    
    G4UIdirectory*             Em6detDir;

    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithAnInteger*      NumOfAbsCmd;

    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsSizYZCmd;

    G4UIcmdWithADoubleAndUnit* AbsXposCmd;

    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldXCmd;
    G4UIcmdWithADoubleAndUnit* WorldYZCmd;

    G4UIcmdWithADoubleAndUnit* MagFieldCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;

};


#endif





