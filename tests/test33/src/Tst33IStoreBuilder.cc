#include "Tst33IStoreBuilder.hh"
#include "G4IStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4GeometryCell.hh"
#include "globals.hh"
#include "Tst33VGeometry.hh"

Tst33IStoreBuilder::Tst33IStoreBuilder()
{}

Tst33IStoreBuilder::~Tst33IStoreBuilder()
{}

G4VIStore *Tst33IStoreBuilder::CreateIStore(Tst33VGeometry *samplegeo) {
  // create an importance store and fill it with the importance
  // per cell values
  const G4VPhysicalVolume &pworld = samplegeo->GetWorldVolume();
  G4IStore *istore = new G4IStore(pworld);
  if (!istore) {
    G4std::G4Exception("Tst33IStoreBuilder::CreateIStore new failed to create G4IStore!");
    }
  // adding GeometryCell for world volume. ReplicaNumer = -1 !
  G4GeometryCell gWorldCell(pworld, -1);
  istore->AddImportanceGeometryCell(1, gWorldCell);
  
  G4int i = 1;
  for (i=1; i <= 19; ++i) {
    G4String volname = samplegeo->GetCellName(i);
    G4double imp = pow(2,i-1);
    if (i==19) {
	imp = pow(2,17);
    }
    const G4VPhysicalVolume *pvol = samplegeo->
      GetPhysicalVolumeByName(volname);
    if (pvol) {
      G4GeometryCell gCell(*pvol, 0);
      istore->AddImportanceGeometryCell(imp, gCell);
    }
  }
  return istore;
}

