#ifndef G4WorldImpMess_hh
#define G4WorldImpMess_hh G4WorldImpMess_hh

#include "G4UImessenger.hh"

class G4ImportanceGeometryConstructor;
class G4UIcmdWithADouble;

class G4WorldImpMess : public G4UImessenger {
public:
  G4WorldImpMess(G4ImportanceGeometryConstructor *);
  void SetNewValue(G4UIcommand * command, G4String newValue);
private:
  G4ImportanceGeometryConstructor *fIGConst;

  G4UIcmdWithADouble *fIbaseCmd;
  G4double fIbase;
  G4bool fIbaseIsSet;
  G4UIcmdWithADouble *fIexpoCmd;
  G4double fIexpo;
  G4bool fIexpoIsSet;

};

#endif
