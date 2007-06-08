#ifndef DetectorSliceDetectorMessenger_h
#define DetectorSliceDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorSliceDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;

class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;


class DetectorSliceDetectorMessenger: public G4UImessenger {

public:

  DetectorSliceDetectorMessenger( DetectorSliceDetectorConstruction* );
  ~DetectorSliceDetectorMessenger();
    
  void SetNewValue( G4UIcommand*, G4String );
    
private:

  DetectorSliceDetectorConstruction* theDetector;

  G4UIdirectory*             theDetectorDir;

  G4UIcmdWithAString*        theTrackerMaterial;
  G4UIcmdWithAString*        theEmAbsorberMaterial;
  G4UIcmdWithAString*        theEmActiveMaterial;
  G4UIcmdWithAString*        theHadAbsorberMaterial;
  G4UIcmdWithAString*        theHadActiveMaterial;
  G4UIcmdWithAString*        theMuonMaterial;

  G4UIcmdWithABool*          theIsEmCalHomogeneous; 
  G4UIcmdWithABool*          theIsHadCalHomogeneous; 
  G4UIcmdWithADouble*        theTrackerLength;
  G4UIcmdWithADouble*        theEmAbsorberTotalLength;
  G4UIcmdWithADouble*        theHadAbsorberTotalLength;
  G4UIcmdWithADouble*        theMuonLength;
  G4UIcmdWithADouble*        theDetectorRadius;
  G4UIcmdWithAnInteger*      theEmActiveLayerNumber;
  G4UIcmdWithADouble*        theEmActiveLayerSize;
  G4UIcmdWithAnInteger*      theHadActiveLayerNumber;
  G4UIcmdWithADouble*        theHadActiveLayerSize;

  G4UIcmdWithoutParameter*   theUpdateCommand;

};

#endif

