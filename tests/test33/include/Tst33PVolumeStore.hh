#ifndef Tst33PVolumeStore_hh
#define Tst33PVolumeStore_hh Tst33PVolumeStore_hh

#include "globals.hh"
#include "g4std/set"
#include "G4GeometryCell.hh"
#include "G4GeometryCellComp.hh"

typedef G4std::set< G4GeometryCell, G4GeometryCellComp > Tst33SetGeometryCell;

class Tst33PVolumeStore {
public:
  Tst33PVolumeStore();
  ~Tst33PVolumeStore();
  
  void AddPVolume(const G4GeometryCell &cell);
  const G4VPhysicalVolume *GetPVolume(const G4String &name) const;
  G4String GetPNames() const;

private:
  Tst33SetGeometryCell fSetGeometryCell;
};



#endif
