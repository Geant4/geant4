
#ifndef TstDrawVox01DetectorMessenger_h
#define TstDrawVox01DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class TstDrawVox01DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;

class TstDrawVox01DetectorMessenger: public G4UImessenger
{
  public:
    TstDrawVox01DetectorMessenger(TstDrawVox01DetectorConstruction * myDC);
    void SetNewValue(G4UIcommand * command,G4String newValues);
  private:
    TstDrawVox01DetectorConstruction * myDetector;
    G4UIdirectory *      mydetDir;
    G4UIcmdWithAString * selDetCmd;
    G4UIcmdWithAString * switchCmd;
    G4UIcmdWithAString * selMatCmd;

};

#endif

