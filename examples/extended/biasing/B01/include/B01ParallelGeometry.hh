#ifndef B01ParallelGeometry_hh
#define B01ParallelGeometry_hh B01ParallelGeometry_hh

#include "B01VGeometry.hh"
#include "B01PVolumeStore.hh"
#include "B01MaterialFactory.hh"


class B01ParallelGeometry : public B01VGeometry {
public:
  B01ParallelGeometry();
  virtual ~B01ParallelGeometry();

  virtual G4VPhysicalVolume &GetWorldVolume() const;
  virtual const G4VPhysicalVolume *GetPhysicalVolumeByName(const G4String& name)
    const ;
  virtual G4String ListPhysNamesAsG4String() const ;
  virtual G4String GetCellName(G4int i); 

private:
  B01ParallelGeometry(const B01ParallelGeometry &);

  void Construct();

  B01ParallelGeometry &operator=(const B01ParallelGeometry &);

  B01MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  B01PVolumeStore fPVolumeStore;

  G4Material *fGalactic;

};



#endif
