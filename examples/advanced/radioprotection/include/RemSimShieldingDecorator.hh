#ifndef RemSimShieldingDecorator_h
#define RemSimShieldingDecorator_h 1

#include "RemSimDecorator.hh"
#include "RemSimVGeometryComponent.hh"
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class RemSimVGeometryComponent;
class RemSimDecorator;

class RemSimShieldingDecorator: public RemSimDecorator
{
public:
  RemSimShieldingDecorator(RemSimVGeometryComponent*);
  ~RemSimShieldingDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 

private:
  void ConstructShielding(G4VPhysicalVolume*);
  G4Box* shielding;
  G4LogicalVolume* shieldingLog;
  G4VPhysicalVolume* shieldingPhys;
  //RemSimVGeometryComponent* component;
};
#endif
