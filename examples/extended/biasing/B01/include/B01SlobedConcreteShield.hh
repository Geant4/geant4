#ifndef B01SlobedConcreteShield_hh
#define B01SlobedConcreteShield_hh B01SlobedConcreteShield_hh

#include "B01VGeometry.hh"
#include "B01PVolumeStore.hh"
#include "B01MaterialFactory.hh"


class B01SlobedConcreteShield : public B01VGeometry {
public:
  B01SlobedConcreteShield();
  ~B01SlobedConcreteShield();

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

  G4Material *fConcrete;
  G4Material *fGalactic;

};



#endif
