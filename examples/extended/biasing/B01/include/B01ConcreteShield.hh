#ifndef B01ConcreteShield_hh
#define B01ConcreteShield_hh B01ConcreteShield_hh

#include "B01VGeometry.hh"
#include "B01PVolumeStore.hh"
#include "B01MaterialFactory.hh"


class B01ConcreteShield : public B01VGeometry {
public:
  B01ConcreteShield();
  ~B01ConcreteShield();

  G4VPhysicalVolume &GetWorldVolume() const;
  const G4VPhysicalVolume *GetPhysicalVolumeByName(const G4String& name)
    const ;
  G4String ListPhysNamesAsG4String() const ;
  
private:
  void Construct();
  B01MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  B01PVolumeStore fPVolumeStore;

  G4Material *fConcrete;
  G4Material *fGalactic;

};



#endif
