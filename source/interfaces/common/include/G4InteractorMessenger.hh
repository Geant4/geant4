// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#ifndef G4InteractorMessenger_h
#define G4InteractorMessenger_h 1

#include "G4UImessenger.hh"

class G4VInteractiveSession;
class G4UIdirectory;
class G4UIcommand;

class G4InteractorMessenger : public G4UImessenger 
{
public:
  G4InteractorMessenger (G4VInteractiveSession* session);
  ~G4InteractorMessenger ();
  void SetNewValue(G4UIcommand* command,G4String newValue);
private:
  G4VInteractiveSession* session;
  G4UIdirectory* interactorDirectory;
  G4UIcommand* addMenu;
  G4UIcommand* addButton;
};

#endif
