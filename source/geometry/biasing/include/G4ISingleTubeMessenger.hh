#ifndef G4ISingleTubeMessenger_hh
#define G4ISingleTubeMessenger_hh G4ISingleTubeMessenger_hh

#include "G4UImessenger.hh"


class G4ITubeFactory;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

class G4ISingleTubeMessenger : public G4UImessenger {
public:
  G4ISingleTubeMessenger(const G4String &,
			 G4ITubeFactory *);
  
  void SetNewValue(G4UIcommand * command, G4String newValue);
 
private:
  
  G4String fCellName;
  G4ITubeFactory *fITubeFactory;
  G4bool fICellCreated;

  G4UIcmdWithADoubleAndUnit *fZminCmd;
  G4double fZmin;
  G4bool fZminIsSet;
  G4UIcmdWithADoubleAndUnit *fZmaxCmd;
  G4double fZmax;
  G4bool fZmaxIsSet;
  G4UIcmdWithADouble *fIbaseCmd;
  G4double fIbase;
  G4bool fIbasesIsSet;
  G4UIcmdWithADouble *fIexpoCmd;
  G4double fIexpo;
  G4bool fIexpoIsSet;

};

#endif
