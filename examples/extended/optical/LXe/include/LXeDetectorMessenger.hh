
#ifndef LXeDetectorMessenger_h
#define LXeDetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class LXeDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithADouble;

class LXeDetectorMessenger: public G4UImessenger
{
public:
  LXeDetectorMessenger(LXeDetectorConstruction*);
  ~LXeDetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  LXeDetectorConstruction*     LXeDetector;
  G4UIdirectory*               detectorDir;
  G4UIdirectory*               volumesDir;
  G4UIcmdWith3VectorAndUnit*   dimensionsCmd;
  G4UIcmdWithADoubleAndUnit*   housingThicknessCmd;
  G4UIcmdWithADoubleAndUnit*   pmtRadiusCmd;
  G4UIcmdWithAnInteger*        nxCmd;
  G4UIcmdWithAnInteger*        nyCmd;
  G4UIcmdWithAnInteger*        nzCmd;
  G4UIcmdWithABool*            sphereCmd;
  G4UIcmdWithADouble*          reflectivityCmd;
  G4UIcmdWithABool*            wlsCmd;
  G4UIcmdWithABool*            lxeCmd;
  G4UIcmdWithAnInteger*        nFibersCmd;
  G4UIcommand*                 updateCmd;
  G4UIcommand*                 defaultsCmd;
  G4UIcmdWithADouble*        MainScintYield;
  G4UIcmdWithADouble*        WLSScintYield;
};

#endif

