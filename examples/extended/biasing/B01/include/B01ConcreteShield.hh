#ifndef B01ConcreteShield_hh
#define B01ConcreteShield_hh B01ConcreteShield_hh

#include "B01VGeometry.hh"
#include "B01PVolumeStore.hh"
#include "B01MaterialFactory.hh"


class B01ConcreteShield : public B01VGeometry {
public:
  B01ConcreteShield();
  virtual ~B01ConcreteShield();

  virtual G4VPhysicalVolume &GetWorldVolume() const;
  virtual const G4VPhysicalVolume *GetPhysicalVolumeByName(const G4String& name)
    const ;
  virtual G4String ListPhysNamesAsG4String() const ;
  
private:
  B01ConcreteShield(const B01ConcreteShield &);
  B01ConcreteShield &operator=(const B01ConcreteShield &);
  void Construct();
  B01MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  B01PVolumeStore fPVolumeStore;

  G4Material *fConcrete;
  G4Material *fGalactic;

};



#endif
