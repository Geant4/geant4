#ifndef RemSimMoonHabitat_h
#define RemSimMoonHabitat_h 1

class RemSimVGeometryComponent;
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class G4VisAttributes;
class RemSimMoonHabitat: public RemSimVGeometryComponent
{
public:
  RemSimMoonHabitat();
  ~RemSimMoonHabitat();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
 
private:
  RemSimMaterial* pMaterial;
};
#endif
