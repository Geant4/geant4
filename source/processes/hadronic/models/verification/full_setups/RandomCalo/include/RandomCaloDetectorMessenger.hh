#ifndef RandomCaloDetectorMessenger_h
#define RandomCaloDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class RandomCaloDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;

class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;


class RandomCaloDetectorMessenger: public G4UImessenger {

public:

  RandomCaloDetectorMessenger( RandomCaloDetectorConstruction* );
  ~RandomCaloDetectorMessenger();
    
  void SetNewValue( G4UIcommand*, G4String );
    
private:

  RandomCaloDetectorConstruction* theDetector;
  G4UIdirectory* theDetectorDir;
  G4UIcmdWithADoubleAndUnit* theFieldCommand;

};

#endif

