#ifndef G4ParameterisedBreast_h
#define G4ParameterisedBreast_h 1

#include "G4VPhysicalVolume.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4HumanPhantomROGeometry;
class G4HumanPhantomSDBreast;

class G4ParameterisedBreast
{
public:
  G4ParameterisedBreast();
  ~G4ParameterisedBreast();
  G4VPhysicalVolume* ConstructBreast(G4VPhysicalVolume*, G4String, G4bool);

private:
  G4HumanPhantomROGeometry* phantomROGeometry;//pointer to ROGeometry
  G4HumanPhantomSDBreast* breastSD;
};
#endif
