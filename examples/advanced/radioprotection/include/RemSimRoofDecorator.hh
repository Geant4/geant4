#ifndef RemSimRoofDecorator_h
#define RemSimRoofDecorator_h 1

#include "RemSimDecorator.hh"
#include "RemSimVGeometryComponent.hh"
#include "globals.hh"
class G4VPhysicalVolume;
class G4Trd;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class RemSimVGeometryComponent;
class RemSimDecorator;
class G4VPhysicalVolume;
class G4VisAttributes;

class RemSimRoofDecorator: public RemSimDecorator
{
public:
  RemSimRoofDecorator(RemSimVGeometryComponent*);
  ~RemSimRoofDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void ChangeMaterial(G4String);
  void PrintDetectorParameters();
  G4VPhysicalVolume* GetShelter(){return 0;};
  void ChangeMother(G4VPhysicalVolume*){;};

private:

  void ConstructRoof(G4VPhysicalVolume*);

  RemSimMaterial* pMaterial; 
  G4String shieldingMaterial;
  G4VisAttributes* roofVisAtt;
  G4Trd* roof;
  G4LogicalVolume* roofLog;
  G4VPhysicalVolume* roofPhys;
};
#endif
