//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....#ifndef 

#ifndef XrayFluoDetectorMessenger_h
#define XrayFluoDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class XrayFluoDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoDetectorMessenger: public G4UImessenger
{
  public:
    XrayFluoDetectorMessenger(XrayFluoDetectorConstruction* );
   ~XrayFluoDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);

 private:
    XrayFluoDetectorConstruction*    Detector;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        SenMaterCmd;
    G4UIcmdWithAString*        SamMaterCmd;  
    G4UIcmdWithADoubleAndUnit* SenThickCmd;
    G4UIcmdWithADoubleAndUnit* SamThickCmd;
    G4UIcmdWithADoubleAndUnit* SenSizeYZCmd;
    G4UIcmdWithADoubleAndUnit* SamSizeYZCmd;
    G4UIcmdWithADoubleAndUnit* SenDistCmd;
    G4UIcmdWithADoubleAndUnit* SenAngleCmd;
    G4UIcmdWithADoubleAndUnit* SenRotCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

#endif
