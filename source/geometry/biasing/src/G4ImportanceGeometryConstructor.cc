#include "G4ImportanceGeometryConstructor.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4UImessenger.hh"
#include "G4Material.hh"
#include "G4IStore.hh"
#include "G4IWorldTubeMessenger.hh"
#include "G4PVPlacement.hh"
#include "G4ImportanceGeometryMessenger.hh"
#include "G4WorldImpMess.hh"

G4ImportanceGeometryConstructor::G4ImportanceGeometryConstructor() :
  fSolidType("none"),
  fSolidMessenger(0),
  fGalactic(new G4Material("Galactic", 
			   1., 
			   1.01*g/mole, 
			   universe_mean_density,
			   kStateGas,
			   2.73*kelvin,
			   3.e-18*pascal)),
  fWorldSolid(0),
  fWorldVolume(0),
  fLogicalWorld(0),
  fIStore(0)
{
  fSolidTypes.insert("tube");
  fGeoMessenger = new G4ImportanceGeometryMessenger(*this);
}
void
G4ImportanceGeometryConstructor::
SetSolidType(const G4String &solidtype) {
  G4ImportanceSolidTypes::iterator it = fSolidTypes.find(solidtype);
  if (it == fSolidTypes.end()){
    Error("Solid not supported");
  }
  fSolidType = *it;
  if (fSolidType=="tube") {
    fSolidMessenger = new G4IWorldTubeMessenger(this);
  } 
}

void
G4ImportanceGeometryConstructor::
ConstructWorldVolume(G4VSolid *worldsolid){
  fWorldSolid = worldsolid;
  
  G4String logname("imp_WorldLog");
  fLogicalWorld = 
    new G4LogicalVolume(fWorldSolid, fGalactic, 
                        logname);

  G4String worldname("ImportanceWorld");
  fWorldVolume = new G4PVPlacement(0, G4ThreeVector(0,0,0), 
				   fLogicalWorld,
				   worldname, 0, false, 0);

  

  fIStore = new G4IStore(*fWorldVolume);
  fIStore->AddImportanceRegion(1, *fWorldVolume, -1); 

  fWorldImpMess = new G4WorldImpMess(this);
  
}

void G4ImportanceGeometryConstructor::
SetWorldImportance(G4double b, G4double e){
  fIStore->ChangeImportance(pow(b,e), *fWorldVolume, -1);
}


G4VIStore *G4ImportanceGeometryConstructor::GetIStore()
{
  return fIStore;
}


G4LogicalVolume *G4ImportanceGeometryConstructor::GetLogicalWorld() {
  return fLogicalWorld;
}




G4VPhysicalVolume *G4ImportanceGeometryConstructor::GetWorldVolume(){
  return fWorldVolume;
}







