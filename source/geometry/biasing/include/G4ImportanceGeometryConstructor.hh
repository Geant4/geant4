#ifndef G4ImportanceGeometryConstructor_hh
#define G4ImportanceGeometryConstructor_hh G4ImportanceGeometryConstructor_hh

#include "globals.hh"
#include "g4std/set"

typedef G4std::set<G4String > G4ImportanceSolidTypes;

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid;
class G4UImessenger;
class G4Material;
class G4VIStore;
class G4ImportanceGeometryMessenger;
class G4WorldImpMess;

class G4ImportanceGeometryConstructor {
public:
  G4ImportanceGeometryConstructor();
  ~G4ImportanceGeometryConstructor(){};
  void SetSolidType(const G4String &solidtype);
  void ConstructWorldVolume(G4VSolid *);
  void SetWorldImportance(G4double b, G4double e);
  G4VIStore *GetIStore();

  G4VPhysicalVolume *GetWorldVolume();
  G4LogicalVolume *GetLogicalWorld();
  
private:

  void Error(const G4String &m){
    G4Exception("Error: G4ImportanceGeometryConstructor: " + m);
  }

  G4ImportanceGeometryMessenger *fGeoMessenger;

  G4String fSolidType;
  G4UImessenger *fSolidMessenger;
  G4Material *fGalactic;
  G4VSolid *fWorldSolid;
  G4VPhysicalVolume *fWorldVolume;
  G4LogicalVolume *fLogicalWorld;
  G4VIStore *fIStore;
  G4ImportanceSolidTypes fSolidTypes;  
  G4WorldImpMess *fWorldImpMess;
};

#endif
