#include "G4WorldImpMess.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4ImportanceGeometryConstructor.hh"

G4WorldImpMess::G4WorldImpMess(G4ImportanceGeometryConstructor *igc) 
{
  fIGConst = igc;

  fIbaseCmd = new G4UIcmdWithADouble("/imp/worldvolume/ibase", this);
  fIbaseIsSet = false;
  fIexpoCmd = new G4UIcmdWithADouble("/imp/worldvolume/iexpo", this);
  fIexpoIsSet = false;

};


void G4WorldImpMess::
SetNewValue(G4UIcommand * command, G4String newValue){
  if (command==fIbaseCmd) {
    fIbase = fIbaseCmd->GetNewDoubleValue(newValue);
    fIbaseIsSet = true;
  }
  if (command==fIexpoCmd) {
    fIexpo = fIexpoCmd->GetNewDoubleValue(newValue);
    fIexpoIsSet = true;
  }
  if (fIbaseIsSet && fIexpoIsSet) {
    fIGConst->SetWorldImportance(fIbase, fIexpo);
  }

}
