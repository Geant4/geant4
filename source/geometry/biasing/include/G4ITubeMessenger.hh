#ifndef G4ITubeMessenger_hh
#define G4ITubeMessenger_hh G4ITubeMessenger_hh

#include "G4UImessenger.hh"
#include "g4std/map"

class G4UIcmdWithAString;
class G4ISingleTubeMessenger;
class G4ITubeFactory;

typedef G4std::map<G4String , G4ISingleTubeMessenger *> G4MapNameTubeMess;

class G4ITubeMessenger : public G4UImessenger {
public:
  G4ITubeMessenger(G4ITubeFactory *);
  void SetNewValue(G4UIcommand * command, G4String newValue);

private:
  void Error(const G4String &m) {
    G4Exception("Error: G4ITubeMessenger" + m);
  }

  G4ITubeFactory *fITubeFactory;

  G4UIcmdWithAString *fCellCreateCmd;
  G4MapNameTubeMess fMapNameTubeMess;

};

#endif
