#ifndef DETECTORMESSENGER_HH
#define DETECTORMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;


class DetectorMessenger : public G4UImessenger {

 public:
   DetectorMessenger(DetectorConstruction*);
   ~DetectorMessenger();

   void SetNewValue(G4UIcommand*, G4String);

 private:
   DetectorConstruction* detectorConstruction;

   G4UIdirectory* targetDirectory;
   G4UIcmdWithAString* frontLayerCmd;
   G4UIcmdWithADoubleAndUnit* layerRadiusCmd;
   G4UIcmdWithADoubleAndUnit* layerThicknessCmd;
   G4UIcmdWithAString* layerMaterialCmd;
   G4UIcmdWithADoubleAndUnit* layerMaxStepSizeCmd;
   G4UIcmdWithADoubleAndUnit* calPositionCmd;
   G4UIcmdWithADoubleAndUnit* calThicknessCmd;
};

#endif // DETECTORMESSENGER_HH
