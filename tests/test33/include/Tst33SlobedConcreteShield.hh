#ifndef Tst33SlobedConcreteShield_hh
#define Tst33SlobedConcreteShield_hh Tst33SlobedConcreteShield_hh

#include "Tst33VGeometry.hh"
#include "Tst33PVolumeStore.hh"
#include "Tst33MaterialFactory.hh"


class Tst33SlobedConcreteShield : public Tst33VGeometry {
public:
  Tst33SlobedConcreteShield();
  virtual ~Tst33SlobedConcreteShield();

  virtual G4VPhysicalVolume &GetWorldVolume() const;
  virtual const G4VPhysicalVolume *GetPhysicalVolumeByName(const G4String& name)
    const ;
  virtual G4String ListPhysNamesAsG4String() const ;
  virtual G4String GetCellName(G4int i); 
  

private:
  Tst33SlobedConcreteShield(const Tst33SlobedConcreteShield &);

  void Construct();

  Tst33SlobedConcreteShield &operator=(const Tst33SlobedConcreteShield &);

  Tst33MaterialFactory fMaterialFactory;
  G4VPhysicalVolume *fWorldVolume;
  Tst33PVolumeStore fPVolumeStore;

  G4Material *fConcrete;
  G4Material *fGalactic;

};



#endif
