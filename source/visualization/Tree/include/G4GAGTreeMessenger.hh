// Satoshi Tanaka 31th May 2001
//
// A messenger for G4GAGTree driver.


#ifndef G4GAGTREEMESSENGER_HH
#define G4GAGTREEMESSENGER_HH

#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4GAGTree;

class G4GAGTreeMessenger: public G4UImessenger {
public:
  G4GAGTreeMessenger(G4GAGTree*);
  virtual ~G4GAGTreeMessenger();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4GAGTree* fpGAGTree;
  G4UIcommand* fpDirectory;
  G4UIcmdWithAnInteger* fpCommandVerbose;
};

#endif
