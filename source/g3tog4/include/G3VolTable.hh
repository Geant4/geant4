#ifndef _G3VOLTABLE_
#define _G3VOLTABLE_ 1

class G4LogicalVolume;
class G4VPhysicalVolume;

class G3VolTable{
private:
  G4LogicalVolume* G3toG4LogicalMother;
  G4VPhysicalVolume* _pv;
  G4int _copy;
  G4VPhysicalVolume* G3toG4PhysicalMother;
  
public:
  G4LogicalVolume* GetLV(const G4String& vname="-");
  G4VPhysicalVolume* GetPV(const G4String& pname="-", G4int pcopy=0);
  void SetMother(G4LogicalVolume* theMother);
  G3VolTable();
  virtual ~G3VolTable();
};
extern G3VolTable G3Vol;
#endif












