
#ifndef Tst05DetectorMessenger_h
#define Tst05DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst05DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst05DetectorMessenger: public G4UImessenger
{
  public:
    Tst05DetectorMessenger(Tst05DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst05DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selDetCmd;
    G4UIcmdWithAString * switchCmd;
    G4UIcmdWithAString * selMatCmd;

};

#endif

