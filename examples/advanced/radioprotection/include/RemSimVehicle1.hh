#ifndef RemSimVehicle1_h
#define RemSimVehicle1_h 1

class RemSimVGeometryComponent;
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class G4VisAttributes;
class RemSimVehicle1: public RemSimVGeometryComponent
{
public:
  RemSimVehicle1();
  ~RemSimVehicle1();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  G4double GetDensity();
  void ChangeMaterial(G4String);
private:
  RemSimMaterial* pMaterial;
  G4VPhysicalVolume* layer1Phys;
};
#endif
