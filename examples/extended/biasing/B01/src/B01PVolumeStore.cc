#include "B01PVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4StringConversion.hh"

B01PVolumeStore::B01PVolumeStore(){}
B01PVolumeStore::~B01PVolumeStore(){}
  
void B01PVolumeStore::AddPVolume(const G4GeometryCell &cell){

  B01SetGeometryCell::iterator it = 
    fSetGeometryCell.find(cell);
  if (it != fSetGeometryCell.end()) {
    G4std::G4cout << "B01PVolumeStore::AddPVolume: cell already stored" 
	   << G4endl;
    return;
  }

  fSetGeometryCell.insert(cell);

    
}

const G4VPhysicalVolume *B01PVolumeStore::
GetPVolume(const G4String &name) const {
  const G4VPhysicalVolume *pvol = 0;
  for (B01SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    if (vol.GetName() == name) {
      pvol =  &vol;
    } 
  }
  if (!pvol) {
    G4std::G4cout << "B01PVolumeStore::GetPVolume: no physical volume named: " 
	   << name << ", found" << G4endl;
  }
  return pvol;
}

G4String B01PVolumeStore::GetPNames() const {
  G4String NameString;
  for (B01SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); ++it) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    G4String cellname(vol.GetName());
    cellname += G4String("_");
    cellname += G4std::str(it->GetReplicaNumber());
    NameString += cellname + G4String("\n");
  }
  return NameString;
}
