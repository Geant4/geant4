//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef exrdmDetectorMessenger_h
#define exrdmDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class exrdmDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class exrdmDetectorMessenger: public G4UImessenger
{
  public:
    exrdmDetectorMessenger(exrdmDetectorConstruction*);
   ~exrdmDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    exrdmDetectorConstruction* myDetector;
    
    G4UIdirectory*             exrdmDir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        TargMatCmd;
    G4UIcmdWithAString*        DetectMatCmd;
    G4UIcmdWithADoubleAndUnit* TargRadiusCmd;
    G4UIcmdWithADoubleAndUnit* DetectThicknessCmd;
    G4UIcmdWithADoubleAndUnit* TargLengthCmd;
    G4UIcmdWithADoubleAndUnit* DetectLengthCmd;
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

