#ifndef RemSimAstronautDecorator_h
#define RemSimAstronautDecorator_h 1

#include "RemSimDecorator.hh"
#include "RemSimVGeometryComponent.hh"
#include "globals.hh"
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class RemSimVGeometryComponent;
class RemSimDecorator;
class G4VPhysicalVolume;
class G4VisAttributes;

class RemSimAstronautDecorator: public RemSimDecorator
{
public:
  RemSimAstronautDecorator(RemSimVGeometryComponent*);
  ~RemSimAstronautDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void ChangeMaterial(G4String);
  void PrintDetectorParameters();

private:
  void ConstructAstronaut(G4VPhysicalVolume*); 
  G4double phantomZ;
  G4double translation;
  G4Box* phantom;
  G4LogicalVolume* phantomLog;
  G4VPhysicalVolume* phantomPhys;
};
#endif
