#ifndef G4ImportanceGeometryMessenger_hh
#define G4ImportanceGeometryMessenger_hh G4ImportanceGeometryMessenger_hh

#include "G4UImessenger.hh"

class G4ImportanceGeometryConstructor;
class G4UIcmdWithAString;

class G4ImportanceGeometryMessenger : public  G4UImessenger
{
public:
  G4ImportanceGeometryMessenger(G4ImportanceGeometryConstructor &igeo);
  void SetNewValue(G4UIcommand * command, G4String newValue);
private:
  
  G4ImportanceGeometryConstructor &fImpGeoConst;
  G4UIcmdWithAString *fSolidTypeCmd;
  
};
#endif
