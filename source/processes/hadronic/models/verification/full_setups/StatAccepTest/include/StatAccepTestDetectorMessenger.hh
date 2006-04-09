#ifndef StatAccepTestDetectorMessenger_h
#define StatAccepTestDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class StatAccepTestDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;

class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;


class StatAccepTestDetectorMessenger: public G4UImessenger {

public:

  StatAccepTestDetectorMessenger( StatAccepTestDetectorConstruction* );
  ~StatAccepTestDetectorMessenger();
    
  void SetNewValue( G4UIcommand*, G4String );
    
private:

  StatAccepTestDetectorConstruction*    theDetector;
  G4UIdirectory*             theDetectorDir;
  G4UIcmdWithADoubleAndUnit* theFieldCommand;
  G4UIcmdWithAString*        theAbsorberMaterial;
  G4UIcmdWithAString*        theActiveMaterial;

  G4UIcmdWithABool*          theIsCalHomogeneous; 
  G4UIcmdWithABool*          theIsUnitInLambda;
  G4UIcmdWithADouble*        theAbsorberTotalLength;
  G4UIcmdWithADouble*        theCalorimeterRadius;
  G4UIcmdWithAnInteger*      theActiveLayerNumber;
  G4UIcmdWithADouble*        theActiveLayerSize;
  G4UIcmdWithAnInteger*      theReadoutLayerNumber;
  // For the specifications of the calorimeter.

  G4UIcmdWithABool*          theIsRadiusUnitInLambda;
  G4UIcmdWithADouble*        theRadiusBinSize;
  G4UIcmdWithAnInteger*      theRadiusBinNumber;
  // For the specifications of the transverse shower analysis.

  G4UIcmdWithoutParameter*   theUpdateCommand;

};

#endif

