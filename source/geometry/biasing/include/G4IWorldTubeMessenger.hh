#ifndef G4IWorldTubeMessanger_hh
#define G4IWorldTubeMessanger_hh G4IWorldTubeMessanger_hh

#include "G4UImessenger.hh"

class G4UIcmdWithADoubleAndUnit;
class G4VSolid;
class G4ITubeFactory;
class G4ImportanceGeometryConstructor;
class G4LogicalVolume;
class G4VIStore;

class G4IWorldTubeMessenger : public  G4UImessenger
{
public:
  G4IWorldTubeMessenger(G4ImportanceGeometryConstructor *);
  void SetNewValue(G4UIcommand * command, G4String newValue);
  void ConstructWorldSolid();
  void ConstructICellFactory(G4LogicalVolume *,
			     G4VIStore *);
private:

  void Error(const G4String& m){
    G4Exception("Error: G4IWorldTubeMessenger: " + m);
  }
  
  G4ImportanceGeometryConstructor *fIGconst;

  G4VSolid *fWorldSolid;
  G4UIcmdWithADoubleAndUnit *fRadiusCmd;
  G4double fRadius;
  G4bool fRadiusIsSet;
  G4UIcmdWithADoubleAndUnit *fHalfHightCmd;
  G4bool fHalfhightIsSet;
  G4double fHalfHight;
  G4ITubeFactory *fITubeFactory;

};
#endif
