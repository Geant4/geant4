#ifndef RemSimShieldingDecorator_h
#define RemSimShieldingDecorator_h 1

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

class RemSimShieldingDecorator: public RemSimDecorator
{
public:
  RemSimShieldingDecorator(RemSimVGeometryComponent*);
  ~RemSimShieldingDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void ChangeMaterial(G4String);
  void PrintDetectorParameters();

private:
  void ConstructShielding(G4VPhysicalVolume*);

  G4double shieldingX;
  G4double shieldingY;
  G4double shieldingZ;
  G4double translation;
  RemSimMaterial* pMaterial; 
  G4String shieldingMaterial;
  G4VisAttributes* shieldingVisAtt;
  G4Box* shielding;
  G4LogicalVolume* shieldingLog;
  G4VPhysicalVolume* shieldingPhys;
};
#endif
