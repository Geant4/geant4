
#ifndef TstVADetectorMessenger_h
#define TstVADetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TstVADetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class TstVADetectorMessenger: public G4UImessenger
{
  public:
    TstVADetectorMessenger(TstVADetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    TstVADetectorConstruction* myDetector;
    G4UIdirectory*             mydetDir;
    G4UIcmdWithAString*        selDetCmd;
    G4UIcmdWithAString*        switchCmd;
    G4UIcmdWithAString*        selMatCmd;

};

#endif

