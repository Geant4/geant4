#ifndef B01PVolumeStore_hh
#define B01PVolumeStore_hh B01PVolumeStore_hh

#include "globals.hh"
#include "g4std/set"
#include "G4GeometryCell.hh"

typedef G4std::set< G4GeometryCell, G4GeometryCellComp > B01SetGeometryCell;

class B01PVolumeStore {
public:
  B01PVolumeStore();
  ~B01PVolumeStore();
  
  void AddPVolume(const G4GeometryCell &cell);
  const G4VPhysicalVolume *GetPVolume(const G4String &name) const;
  G4String GetPNames() const;

private:
  B01SetGeometryCell fSetGeometryCell;
};

#endif
