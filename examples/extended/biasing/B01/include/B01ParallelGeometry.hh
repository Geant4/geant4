#ifndef B01ParallelGeometry_hh
#define B01ParallelGeometry_hh B01ParallelGeometry_hh

#include "B01VGeometry.hh"
#include "B01PVolumeStore.hh"
#include "B01MaterialFactory.hh"


class B01ParallelGeometry : public B01VGeometry {
public:
  B01ParallelGeometry();
  ~B01ParallelGeometry();

  G4VPhysicalVolume &GetWorldVolume() const;
  const G4VPhysicalVolume *GetPhysicalVolumeByName(const G4String& name)
    const ;
  G4String ListPhysNamesAsG4String() const ;
  G4String GetCellName(G4int i); 

private:
  void Construct();
  B01MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  B01PVolumeStore fPVolumeStore;

  G4Material *fGalactic;

};



#endif
