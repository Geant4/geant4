
#ifndef Tst01DetectorMessenger_h
#define Tst01DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst01DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst01DetectorMessenger: public G4UImessenger
{
  public:
    Tst01DetectorMessenger(Tst01DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst01DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selDetCmd;
    G4UIcmdWithAString * switchCmd;
    G4UIcmdWithAString * selMatCmd;

};

#endif

