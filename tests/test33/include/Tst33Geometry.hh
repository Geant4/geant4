#ifndef Tst33VGeometry_hh
#define Tst33VGeometry_hh Tst33VGeometry_hh

#include "globals.hh"

class G4VPhysicalVolume;

class Tst33VGeometry {
public:
  Tst33VGeometry();
  virtual ~Tst33VGeometry();

  const G4VPhysicalVolume &GetWorldVolume() const;
  const G4VPhysicalVolume &GetPhysicalVolumeByName(const G4String& name)
    const ;
  G4String ListPhysNamesAsG4String() const ;
  
}


#endif
