#ifndef RemSimAstronautDecorator_h
#define RemSimAstronautDecorator_h 1

#include "RemSimDecorator.hh"
#include "RemSimVGeometryComponent.hh"
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class RemSimVGeometryComponent;
class RemSimDecorator;

class RemSimAstronautDecorator: public RemSimDecorator
{
public:
  RemSimAstronautDecorator(RemSimVGeometryComponent*);
  ~RemSimAstronautDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 

private:
  void ConstructAstronaut(G4VPhysicalVolume*);
  G4Box* astronaut;
  G4LogicalVolume* astronautLog;
  G4VPhysicalVolume* astronautPhys;
  //RemSimVGeometryComponent* component;
};
#endif
