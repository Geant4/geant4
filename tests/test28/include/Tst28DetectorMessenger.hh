#ifndef Tst28DetectorMessenger_h
#define Tst28DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst28DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class Tst28DetectorMessenger: public G4UImessenger
{
  public:
    Tst28DetectorMessenger(Tst28DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    Tst28DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selMatCmd;

};

#endif

