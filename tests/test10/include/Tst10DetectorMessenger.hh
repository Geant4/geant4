
#ifndef Tst10DetectorMessenger_h
#define Tst10DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst10DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst10DetectorMessenger: public G4UImessenger
{
  public:
    Tst10DetectorMessenger(Tst10DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst10DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selDetCmd;
};

#endif

