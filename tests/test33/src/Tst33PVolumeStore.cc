#include "Tst33PVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4StringConversion.hh"

Tst33PVolumeStore::Tst33PVolumeStore(){}
Tst33PVolumeStore::~Tst33PVolumeStore(){}
  
void Tst33PVolumeStore::AddPVolume(const G4GeometryCell &cell){

  Tst33SetGeometryCell::iterator it = 
    fSetGeometryCell.find(cell);
  if (it != fSetGeometryCell.end()) {
    G4std::G4cout << "Tst33PVolumeStore::AddPVolume: cell already stored" 
	   << G4endl;
    return;
  }

  fSetGeometryCell.insert(cell);

    
}

const G4VPhysicalVolume *Tst33PVolumeStore::
GetPVolume(const G4String &name) const {
  const G4VPhysicalVolume *pvol = 0;
  for (Tst33SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    if (vol.GetName() == name) {
      pvol =  &vol;
    } 
  }
  if (!pvol) {
    G4std::G4cout << "Tst33PVolumeStore::GetPVolume: no physical volume named: " 
	   << name << ", found" << G4endl;
  }
  return pvol;
}

G4String Tst33PVolumeStore::GetPNames() const {
  G4String NameString;
  for (Tst33SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    G4String cellname(vol.GetName());
    cellname += G4String("_");
    cellname += G4std::str(it->GetReplicaNumber());
    NameString += cellname + G4String("\n");
  }
  return NameString;
}
