#ifndef Tst33VGeometry_hh
#define Tst33VGeometry_hh Tst33VGeometry_hh

#include "globals.hh"

class G4VPhysicalVolume;

class Tst33VGeometry {
public:
  Tst33VGeometry();
  virtual ~Tst33VGeometry();

  virtual G4VPhysicalVolume &GetWorldVolume() const = 0;

  virtual const G4VPhysicalVolume *
  GetPhysicalVolumeByName(const G4String& name) const = 0;
  virtual G4String GetCellName(G4int i) = 0;

  virtual G4String ListPhysNamesAsG4String() const = 0;
  
};


#endif
