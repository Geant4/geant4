#include "B01PVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Pstring.hh"

B01PVolumeStore::B01PVolumeStore(){}
B01PVolumeStore::~B01PVolumeStore(){}
  
void B01PVolumeStore::AddPVolume(const G4GeometryCell &cell){

  B01SetGeometryCell::iterator it = 
    fSetGeometryCell.find(cell);
  if (it != fSetGeometryCell.end()) {
    G4cout << "B01PVolumeStore::AddPVolume: cell already stored" 
	   << G4endl;
    return;
  }

  fSetGeometryCell.insert(cell);

    
}

const G4VPhysicalVolume *B01PVolumeStore::
GetPVolume(const G4String &name) const {
  for (B01SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); it++) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    if (vol.GetName() == name) {
      return &vol;
    } 
  }
  G4cout << "B01PVolumeStore::GetPVolume: no physical volume named: " 
	 << name << ", found" << G4endl;
  return 0;
}

G4String B01PVolumeStore::GetPNames() const {
  G4String NameString;
  for (B01SetGeometryCell::const_iterator it = fSetGeometryCell.begin();
       it != fSetGeometryCell.end(); it++) {
    const G4VPhysicalVolume &vol = it->GetPhysicalVolume();
    G4String cellname;
    cellname = vol.GetName() + "_" + str(it->GetReplicaNumber());
    NameString += cellname + "\n";
  }
  return NameString;
}
