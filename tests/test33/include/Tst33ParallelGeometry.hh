#ifndef Tst33ParallelGeometry_hh
#define Tst33ParallelGeometry_hh Tst33ParallelGeometry_hh

#include "Tst33VGeometry.hh"
#include "Tst33PVolumeStore.hh"
#include "Tst33MaterialFactory.hh"


class Tst33ParallelGeometry : public Tst33VGeometry {
public:
  Tst33ParallelGeometry();
  virtual ~Tst33ParallelGeometry();

  virtual G4VPhysicalVolume &GetWorldVolume() const;
  virtual const G4VPhysicalVolume *GetPhysicalVolumeByName(const G4String& name)
    const ;
  virtual G4String ListPhysNamesAsG4String() const ;
  virtual G4String GetCellName(G4int i); 

private:
  Tst33ParallelGeometry(const Tst33ParallelGeometry &);

  void Construct();

  Tst33ParallelGeometry &operator=(const Tst33ParallelGeometry &);

  Tst33MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  Tst33PVolumeStore fPVolumeStore;

  G4Material *fGalactic;

};



#endif
