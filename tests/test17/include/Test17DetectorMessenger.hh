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

#ifndef Test17DetectorMessenger_h
#define Test17DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Test17DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17DetectorMessenger: public G4UImessenger
{
public: // Without description

    Test17DetectorMessenger(Test17DetectorConstruction* );
   ~Test17DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Test17DetectorConstruction*   Test17Detector;
    
    G4UIdirectory*             Test17detDir;

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





