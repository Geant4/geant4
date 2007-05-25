#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
  public:
  
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
    DetectorConstruction*   Detector;
    
    G4UIdirectory*             testemDir;
    G4UIdirectory*             detDir;    
    G4UIcmdWithAString*        Mater1Cmd;
    G4UIcmdWithAString*        Mater2Cmd;
    G4UIcmdWithAString*        Mater3Cmd;
    G4UIcmdWithADoubleAndUnit* SizeX1Cmd;
    G4UIcmdWithADoubleAndUnit* SizeX2Cmd;
    G4UIcmdWithADoubleAndUnit* SizeX3Cmd;
    G4UIcmdWithAnInteger*      NbLayers1Cmd;        
    G4UIcmdWithAnInteger*      NbLayers2Cmd;        
    G4UIcmdWithAnInteger*      NbLayers3Cmd;        
            
    G4UIcmdWithoutParameter*   UpdateCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

