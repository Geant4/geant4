#ifndef RemSimShelterSPEDecorator_h
#define RemSimShelterSPEDecorator_h 1

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

class RemSimShelterSPEDecorator: public RemSimDecorator
{
public:
  RemSimShelterSPEDecorator(RemSimVGeometryComponent*);
  ~RemSimShelterSPEDecorator();

  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void ChangeMaterial(G4String);
  void PrintDetectorParameters();
  G4VPhysicalVolume* GetShelter(){return 0;};
  void ChangeMother(G4VPhysicalVolume*){;};

private:
  void ConstructShelterSPE(G4VPhysicalVolume*);

  G4double shelterSPEX;
  G4double shelterSPEY;
  G4double shelterSPEZ;
  G4double translation;
  RemSimMaterial* pMaterial; 
  G4String shelterSPEMaterial;
  G4VisAttributes* shelterSPEVisAtt;
  G4Box* shelterSPE;
  G4LogicalVolume* shelterSPELog;
  G4VPhysicalVolume* shelterSPEPhys;
};
#endif
