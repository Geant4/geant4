#include "G4ISingleTubeMessenger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4ITubeFactory.hh"

G4ISingleTubeMessenger::G4ISingleTubeMessenger(const G4String &cellname,
					       G4ITubeFactory *tfact) {

  fCellName = cellname;
  fITubeFactory = tfact;
  fZminIsSet = false;
  fZmaxIsSet = false;
  fIbasesIsSet = false;
  fIexpoIsSet = false;
  fICellCreated = false;
  
  G4String path("/imp/cell/" + fCellName + "/");

  fZminCmd = new G4UIcmdWithADoubleAndUnit(G4String(path + "zmin"),this);
  fZmaxCmd = new G4UIcmdWithADoubleAndUnit(G4String(path + "zmax"),this);
  fIbaseCmd = new G4UIcmdWithADouble(G4String(path + "ibase"),this);
  fIexpoCmd = new G4UIcmdWithADouble(G4String(path + "iexpo"),this);

}

void G4ISingleTubeMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){

  if (command==fZminCmd) {
    fZmin = fZminCmd->GetNewDoubleValue(newValue);
    fZminIsSet = true;
  }
  if (command==fZmaxCmd) {
    fZmax = fZmaxCmd->GetNewDoubleValue(newValue);
    fZmaxIsSet = true;
  }
  if (command==fIbaseCmd) {
    fIbase = fIbaseCmd->GetNewDoubleValue(newValue);
    fIbasesIsSet = true;
  }
  if (command==fIexpoCmd) {
    fIexpo = fIexpoCmd->GetNewDoubleValue(newValue);
    fIexpoIsSet = true;
  }
  
  if (fZminIsSet && fZmaxIsSet && fIbasesIsSet && fIexpoIsSet) {
    if (!fICellCreated) {
      fITubeFactory->AddCell(fCellName, fZmin, fZmax);
      fICellCreated = true;
    }
    fITubeFactory->SetImportance(fCellName, fIbase, fIexpo);
  }

}



