#include "G4IWorldTubeMessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Tubs.hh"
#include "G4ImportanceGeometryConstructor.hh"
#include "G4LogicalVolume.hh"
#include "G4VIStore.hh"
#include "G4ITubeFactory.hh"

G4IWorldTubeMessenger::
G4IWorldTubeMessenger(G4ImportanceGeometryConstructor *igc):
  fIGconst(igc){
  fRadiusCmd = 
    new G4UIcmdWithADoubleAndUnit("/imp/worldvolume/tube/radius", this);
  fRadiusIsSet =false;
  fHalfHightCmd = 
    new G4UIcmdWithADoubleAndUnit("/imp/worldcolume/tube/halfwidth", this);
  fHalfhightIsSet = false;
  fITubeFactory = 0;
}

void G4IWorldTubeMessenger::
SetNewValue(G4UIcommand * command, G4String newValue){
  
  if (command==fRadiusCmd) {
    fRadius = fRadiusCmd->GetNewDoubleValue(newValue);
    fRadiusIsSet = true;
  }
  if (command==fHalfHightCmd) {
    fHalfHight = fHalfHightCmd->GetNewDoubleValue(newValue);
    fHalfhightIsSet = true;
  }
  
  if (fRadiusIsSet && fHalfhightIsSet) {
    ConstructWorldSolid();
    fIGconst->ConstructWorldVolume(fWorldSolid);
    ConstructICellFactory(fIGconst->GetLogicalWorld(),
			  fIGconst->GetIStore());
  }
  
}

void G4IWorldTubeMessenger::ConstructWorldSolid() {
  if (!fRadiusIsSet || !fHalfhightIsSet) {
    Error("!fRadiusIsSet || !fHalfhightIsSet");
  }
  if (!fWorldSolid) {
    G4String tubename("imp_WorldTube");
    fWorldSolid = new G4Tubs(tubename,
			     0.,
			     fRadius,
			     fHalfHight,
			     0.*deg,
			     360.*deg);
    G4cout << "G4IWorldTubeMessenger:: constructed World solid, with" 
	   << G4endl
	   <<"   radius    = " << fRadius << G4endl
	   <<"   halfhight = " 
	   << fHalfHight  << G4endl;
  }
}

void G4IWorldTubeMessenger::ConstructICellFactory(G4LogicalVolume *lv,
						  G4VIStore *is){
  if (!fRadiusIsSet || !fHalfhightIsSet) {
    Error("!fRadiusIsSet || !fHalfhightIsSet");
  }
  if (!fITubeFactory) {
    fITubeFactory = new G4ITubeFactory(lv,
				       is,
				       fRadius, 
				       fHalfHight);
  }
}

