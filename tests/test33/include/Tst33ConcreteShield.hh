#ifndef Tst33ConcreteShield_hh
#define Tst33ConcreteShield_hh Tst33ConcreteShield_hh

#include "Tst33VGeometry.hh"
#include "Tst33PVolumeStore.hh"
#include "Tst33MaterialFactory.hh"


class Tst33ConcreteShield : public Tst33VGeometry {
public:
  Tst33ConcreteShield();
  virtual ~Tst33ConcreteShield();

  virtual G4VPhysicalVolume &GetWorldVolume() const;
  virtual const G4VPhysicalVolume *GetPhysicalVolumeByName(const G4String& name)
    const ;
  virtual G4String ListPhysNamesAsG4String() const ;
  virtual G4String GetCellName(G4int i); 
  
private:
  Tst33ConcreteShield(const Tst33ConcreteShield &);
  Tst33ConcreteShield &operator=(const Tst33ConcreteShield &);
  void Construct();
  Tst33MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  Tst33PVolumeStore fPVolumeStore;

  G4Material *fConcrete;
  G4Material *fGalactic;

};



#endif
