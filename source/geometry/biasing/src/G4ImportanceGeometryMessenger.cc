#include "G4ImportanceGeometryMessenger.hh"
#include "G4ImportanceGeometryConstructor.hh"
#include "G4UIcmdWithAString.hh"

G4ImportanceGeometryMessenger::
G4ImportanceGeometryMessenger(G4ImportanceGeometryConstructor &igeo):
  fImpGeoConst(igeo)
{
  fSolidTypeCmd = new G4UIcmdWithAString("/imp/worldvolume/solid",
					 this);
}

void G4ImportanceGeometryMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){
  
  if (command==fSolidTypeCmd) {
    fImpGeoConst.SetSolidType(newValue);
  }
  
}
