#ifndef RemSimVehicle2_h
#define RemSimVehicle2_h 1

#include "globals.hh"

class RemSimVGeometryComponent;
class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class G4VisAttributes;
class RemSimVehicle2: public RemSimVGeometryComponent
{
public:
  RemSimVehicle2();
  ~RemSimVehicle2();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  G4double GetDensity(){return 0;};
  void ChangeMaterial(G4String){ G4cout << 
                                 "This command is not available for Vehicle2"
					<<G4endl; }
private:
  RemSimMaterial* pMaterial;
  G4Tubs* layer1;
  G4Tubs* layer2;
  G4Tubs* layer3;
  G4Tubs* layer4;
  G4Tubs* layer5;
  G4LogicalVolume* layer1Log;
  G4LogicalVolume* layer2Log;
  G4LogicalVolume* layer3Log;
  G4LogicalVolume* layer4Log;
  G4LogicalVolume* layer5Log;
  G4VPhysicalVolume* layer1Phys;
  G4VPhysicalVolume* layer2Phys;
  G4VPhysicalVolume* layer3Phys;
  G4VPhysicalVolume* layer4Phys;
  G4VPhysicalVolume* layer5Phys;
  G4VisAttributes* layer1VisAtt;
  G4VisAttributes* layer2VisAtt;
  G4VisAttributes* layer3VisAtt;
  G4VisAttributes* layer4VisAtt;
  G4VisAttributes* layer5VisAtt;
};
#endif
