#ifndef MLColourMessenger_h
#define MLColourMessenger_h 1
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UImessenger.hh"

class MLColour;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
////////////////////////////////////////////////////////////////////////////////
//
class MLColourMessenger: public G4UImessenger
{
public:
  MLColourMessenger(MLColour* );
  ~MLColourMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:

  MLColour                  *colourManager;

  G4UIdirectory             *ColourDir;
  G4UIcmdWithoutParameter   *ListCmd;
  G4UIcommand               *AddCmd;
};
////////////////////////////////////////////////////////////////////////////////
#endif
