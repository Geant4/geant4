#ifndef G4ITubeFactory_hh 
#define G4ITubeFactory_hh  G4ITubeFactory_hh 

#include "globals.hh"
#include "g4std/map"

class G4ITubeMessenger;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VIStore;
class G4Material;
typedef G4std::map<G4String, G4Tubs *>  G4MapNameTube;
typedef G4std::map<G4String, G4LogicalVolume *>  G4MapNameLogic;
typedef G4std::map<G4String, G4VPhysicalVolume *>  G4MapNamePhysical;


class G4ITubeFactory {
public:
  G4ITubeFactory(G4LogicalVolume *wl, G4VIStore *is,
		 G4double Radius, G4double HalfHight);
  ~G4ITubeFactory();
  void AddCell(const G4String &celname, G4double zmin, G4double zmax);
  void SetImportance(const G4String &celname, G4double b, G4double e);

private:

  void Error(const G4String &m) {
    G4Exception("Error: G4ITubeFactory: " + m);
  }

  G4LogicalVolume *fWorldLogic;
  G4VIStore *fIStore;
  G4double fRadius;
  G4double fHalfHight;
  G4ITubeMessenger *fITubeMessenger;
  G4Material *fGalactic;

  G4MapNameTube fMapNameTube;
  G4MapNameLogic fMapNameLogic;
  G4MapNamePhysical fMapNamePhysical;


}; 

#endif
