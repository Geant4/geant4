#ifndef B01VGeometry_hh
#define B01VGeometry_hh B01VGeometry_hh

#include "globals.hh"

class G4VPhysicalVolume;

class B01VGeometry {
public:
  virtual ~B01VGeometry(){}

  const G4VPhysicalVolume &GetWorldVolume() const;
  const G4VPhysicalVolume &GetPhysicalVolumeByName(const G4String& name)
    const ;
  G4String ListPhysNamesAsG4String() const ;
  
}


#endif
