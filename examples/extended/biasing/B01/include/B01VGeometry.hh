#ifndef B01VGeometry_hh
#define B01VGeometry_hh B01VGeometry_hh

#include "globals.hh"

class G4VPhysicalVolume;

class B01VGeometry {
public:
  B01VGeometry();
  virtual ~B01VGeometry();

  virtual G4VPhysicalVolume &GetWorldVolume() const = 0;

  virtual const G4VPhysicalVolume *
  GetPhysicalVolumeByName(const G4String& name) const = 0;

  virtual G4String ListPhysNamesAsG4String() const = 0;
  
};


#endif
