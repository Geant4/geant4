#ifndef MyDetectorMessenger_h
#define MyDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MyDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;


class MyDetectorMessenger: public G4UImessenger {

public:

  MyDetectorMessenger(MyDetectorConstruction*);
  ~MyDetectorMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  MyDetectorConstruction*    theDetector;
  G4UIdirectory*             theDetectorDir;
  G4UIcmdWithADoubleAndUnit* theFieldCommand;
  G4UIcmdWithAString*        theAbsorberMaterial;
  G4UIcmdWithAString*        theActiveMaterial;

};

#endif

