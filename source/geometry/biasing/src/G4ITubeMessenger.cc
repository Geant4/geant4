#include "G4ITubeMessenger.hh"
#include "G4ITubeFactory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ISingleTubeMessenger.hh"


G4ITubeMessenger::G4ITubeMessenger(G4ITubeFactory *itfact) :
  fITubeFactory(itfact)
{  
  fCellCreateCmd = new G4UIcmdWithAString("/imp/cell/create", this);
}

void G4ITubeMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){
  if (command==fCellCreateCmd) {
    G4MapNameTubeMess::iterator it = fMapNameTubeMess.find(newValue);
    if (it!=fMapNameTubeMess.end()) {
      Error("cell: " + newValue + ", already exists!");
    }
    fMapNameTubeMess[newValue] = new G4ISingleTubeMessenger(newValue, 
							fITubeFactory);
  }
}

